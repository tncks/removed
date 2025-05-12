
package modi;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

import msutil.IonGraph;
import msutil.PGraph;
import msutil.Scoring;

import processedDB.ProtDatabase;

public class PtmMap {

	private final int[][] ptmFreqTable;
	private final ArrayList<PSM> history= new ArrayList<>();
	
	private final int lowestCol = (int)Math.round( Constants.minModifiedMass );
	private final int highestCol = (int)Math.round( Constants.maxModifiedMass );
	private final int columnSize= highestCol - lowestCol + 10;
	private final int siteColumn= 23;
	
	public PtmMap(){			
		ptmFreqTable= new int[columnSize][siteColumn];		
	}
	public void removeLast( ScanCap scan, double cut ){ 
		int last = history.size();
		if( last < 1 ) return;
		
		PSM p = history.get(last-1);	
		if( scan == p.spec && cut >= p.prob ){
			ptmFreqTable[p.delta-lowestCol][p.getModRes()+2]--;
			if( p.site == 0 ) ptmFreqTable[p.delta-lowestCol][0]--;			
			history.remove(last-1);
		}
	}
	
	public void hitPtmFreqTable( ScanCap spec, Peptide peptide, PTM ptm, PGraph pg ){
	
		double mass= ptm.getMassDifference(); 
		AminoAcid aa= ptm.getResidue();		
		String residue = ( ptm.getPTMPosition() == PTMPosition.ANY_N_TERM )? "nterm" : aa.toString();
		
		int delta = (int)Math.round(mass);
		if( delta < lowestCol || highestCol < delta ) return;
	
		ptmFreqTable[delta-lowestCol][aa.getIndex()+2]++;
		if( residue.compareToIgnoreCase("nterm") == 0 )
			ptmFreqTable[delta-lowestCol][0]++;
		
		double[] ptms = new double[peptide.size()];
		ptms[ptm.getID()] = delta;
		IonGraph iG = Scoring.PeptideSpectrumMatch(peptide.toString(), ptms, pg );
		
		history.add( new PSM(spec, peptide, ptm.getID(), delta, iG.getProb()) );
	}
	
	public void getPtmMap(String source, ProtDatabase ProtDB) throws IOException{		
	
		int i;
		System.out.println("Print PTM Map!");
		
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(source+".modmap.csv")));
		out.println("[PTM Frequency Matrix]");
		out.println("Delta|Nterm|  A  |  C  |  D  |  E  |  F  |  G  |  H  |  I  |  K  |  L  |  M  |  N  |  P  |  Q  |  R  |  S  |  T  |  V  |  W  |  Y  |");
		for(i=0; i<columnSize; i++){
			int hit=0;
			for(int j=0; j<siteColumn; j++) {
				if( j==1 || j==22 ) continue;
				hit += ptmFreqTable[i][j];
			}
			if( hit == 0 )
				continue;
			
			out.print( String.format("%5d", (i+lowestCol)) + "|" );
			for(int k=0; k<siteColumn; k++) {
				if( k==1 || k==22 ) continue;
				out.print( String.format("%5d", ptmFreqTable[i][k]) + "|" );
			}
			out.println();
		}


		Collections.sort(history);
		out.println();
		out.println("[Peptide-Spectrum Matches]");
		out.println("Spectrum_name,"+
					"Observed_M/Z," +
					"CS," 			+
					"CalcMW," 		+
					"DeltaM," 		+
					"Probability," 	+
					"Protein," 		+	
					"Start_site," 	+		
					"End_site," 	+		
					"Peptide,"		);
		
		
		for(int k=0; k<history.size(); k++)
			out.println( history.get(k).toIdentification(ProtDB) );
		
		out.close();
		System.out.println("Done PTM Map!");
	}	
	public void getMODMapToPRIX( PrintWriter out ) {
		
		int i;
		System.out.print("Generating MODMap.....  ");
		
		int candidate = 0;
		for(i=0; i<columnSize; i++){ // delta
			int hit=0;
			for(int j=0; j<siteColumn; j++) { // residue
				if( j==1 || j==22 ) continue;
				hit += ptmFreqTable[i][j];
			}
			if( hit != 0 ) candidate++;
		}
		out.print( candidate+ "|" );
		for(i=0; i<columnSize; i++){
			int hit=0;
			for(int j=0; j<siteColumn; j++) {
				if( j==1 || j==22 ) continue;
				hit += ptmFreqTable[i][j];
			}
			if( hit == 0 ) continue;
			
			out.print( (i+lowestCol)+ "|" );
			for(int k=0; k<siteColumn; k++) {
				if( k==1 || k==22 ) continue;
				out.print( ptmFreqTable[i][k]+ "|" );
			}
		}
		System.out.println( candidate+" modifications detected" );	
	}
	
	private static class PSM implements Comparable<PSM> {
		final ScanCap spec;
		final Peptide pept;
		final int site;
		final int delta;
		final double prob;
		public PSM( ScanCap s, Peptide p, int t, int dt, double score ){
			spec = s;
			pept = p;
			site = t;
			delta = dt;
			prob = score;
		}
		private int getModRes() { return pept.get(site).getIndex(); }
		
		public String toIdentification(ProtDatabase protDB){
			StringBuffer des = new StringBuffer();
			des.append( String.format("%s,%.4f,%d,", spec.getTitle(), spec.getObservedMW(), spec.getCharge()) );
			
			double mw = pept.getMolecularMass()+delta;
			des.append( String.format("%.4f,%.4f,%.4f,", mw, (spec.getObservedMW()-mw), prob) );
			
			des.append( protDB.getProteinIdentity(pept)+"," );
			int start = protDB.getRealPositionInProtein(pept.getSrcProteinInfo().getFirst().getStartPos());
			des.append( String.format("%d,%d,", start, start+pept.size()-1) );
			
			StringBuffer mp = new StringBuffer( protDB.getCuratedPeptide(pept) );
			mp.insert(site+3, delta);
			if( delta > 0 ) mp.insert(site+3, "+");
			des.append(mp);
			
			return des.toString();
		}

		public int compareTo(PSM p){
			
			if( this.getModRes() >  p.getModRes() ) return 1;
			else if( this.getModRes() <  p.getModRes() ) return -1;
			
			if( this.delta > p.delta ) return 1;
			else if( this.delta < p.delta ) return -1;
			else return 0;
		}
	}
	
}
