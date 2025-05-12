package scaniter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import modi.Constants;
import modi.Peak;
import modi.Spectrum;
import msutil.MSMass;

public class MSMScan {
	
	static int minPeaksCount = 4;
	static double minMW = 8 * MSMass.getMinAAMass() + 18;
	
	private String 		title;
	private int 		specIndex;
	private int 		scanNo;
	private double 		pmz;
	private double 		neutralMW;
	private int 		charge;
	private long 		offset;
	private Spectrum 	peaklist;
	private static double tolerance= Constants.massToleranceForDenovo;
	
	private double 	fragmentTol = 0;
	private double 	precursorTolerance = 0;
	private double 	precursorAccuracy= 0;
	private double 	gapTolerance = 0;
	private double 	nonModifiedDelta = 0;
	private int		maxNoOfC13 = 0;
	
	public MSMScan(int index, double pmz, int charge){
		this.title 		= "";
		this.specIndex	= index;
		this.scanNo		= 0;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public MSMScan(String title, int index, int sn, double pmz, int charge){		
		this.title 		= title;
		this.specIndex  = index;
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public String 	getTitle(){ return title; }
	public int 		getSpecIndex(){ return specIndex; }
	public int 		getScanNumber(){ return scanNo; }
	public double 	getObservedMW(){ return neutralMW; }
	public double 	getPMZ(){ return pmz; }
	public int 		getCharge(){ return charge; }
	public int 		getPeaksCount(){ return peaklist.size(); }
	public long 	getOffset(){ return offset; }	
	public String 	getHeader(){ return String.format("%d\t%.4f\t%d\t%d\t%s",
													specIndex, neutralMW, charge, scanNo, title); }
	
	public Spectrum getSpectrum() { 
		Constants.precursorTolerance= precursorTolerance;
		Constants.precursorAccuracy	= precursorAccuracy;
		Constants.gapTolerance 		= gapTolerance;
		Constants.gapAccuracy 		= precursorAccuracy + 2*Constants.fragmentTolerance;
		Constants.nonModifiedDelta 	= nonModifiedDelta;
		Constants.maxNoOfC13 		= maxNoOfC13;
		return peaklist; 
	}
	
	public double 	getPrecursorTolerance() { return precursorTolerance; }
	public double 	getPrecursorAccuracy() { return precursorAccuracy; }
	public double	getGapTolerance() { return gapTolerance; }
	public double 	getGapAccuracy() { return precursorAccuracy + 2*Constants.fragmentTolerance; }
	public double 	getNonModifiedDelta() { return nonModifiedDelta; }
	public int 		getMaxNoOfC13() { return maxNoOfC13; }
	
	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {
		
		if( neutralMW < minMW || Constants.maxPeptideMass < neutralMW ) return false; 
				
		if( Constants.reporterMassOfIsobaricTag != null ) removeReporterIons(rawPL, Constants.reporterMassOfIsobaricTag);
		
		if( Constants.rangeForIsotopeIncrement != 0 ) maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );
		else maxNoOfC13 = Constants.maxNoOfC13;
		
		if( Constants.PPMTolerance != 0 ) precursorAccuracy = Constants.PPMtoDalton( neutralMW, Constants.PPMTolerance );
		else precursorAccuracy = Constants.precursorAccuracy;
		
		precursorTolerance = precursorAccuracy + maxNoOfC13*Constants.IsotopeSpace;
		
		int index = 0;
		Spectrum spectrum = new Spectrum( this.pmz, this.charge, this.title );
		
		double basePeakIntensity=0, TIC=0;
		double tarMass=0, tarInten=0;
		for( RawPeak rp : rawPL ) {		
			double mass = rp.mz;
			double intensity = rp.it;
			if( intensity <= 0 || mass <= 0 ) continue;
			if( mass > neutralMW ) continue;
			
			if( ( mass - tarMass ) < tolerance ){
				double sum = tarInten + intensity;
				tarMass = tarMass*(tarInten/sum)+ mass*(intensity/sum);
				tarInten += intensity;
				spectrum.get(index-1).set(tarMass, tarInten);
			}
			else{
				spectrum.add( new Peak(index++, mass, intensity) );
				tarMass = mass;
				tarInten = intensity; 
			}
			TIC += intensity;
			if( tarInten > basePeakIntensity )
				basePeakIntensity= tarInten;
		}	
		spectrum.setExtraInformation( basePeakIntensity, TIC );
		
		gapTolerance = Constants.fragmentTolerance*2;
		nonModifiedDelta = (precursorTolerance < Constants.massToleranceForDenovo)? precursorTolerance : Constants.massToleranceForDenovo;
				
		if( precursorTolerance > gapTolerance ) gapTolerance += precursorTolerance;
		
		if( spectrum.size() < minPeaksCount ) peaklist = null; 
		else peaklist = spectrum;
		
		return (peaklist!=null);
	}
	
	private void removeReporterIons( ArrayList<RawPeak> rawPL, double[] removedMasses ){
	
		ArrayList<RawPeak> reporters = new ArrayList<RawPeak>();
		for(int i=1; i<removedMasses.length; i++)
			reporters.add( new RawPeak(removedMasses[i], Constants.fragmentTolerance) );
		
		reporters.add( new RawPeak(removedMasses[0] + Constants.Proton, Constants.fragmentTolerance) );

		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - removedMasses[0] + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawPeak(compItraqTag, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(secondIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(thirdIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(forthIso, Constants.fragmentTolerance) );

			for(int i=1; i<=maxNoOfC13; i++){
				reporters.add( new RawPeak(compItraqTag-i*Constants.IsotopeSpace/fragCS, Constants.fragmentTolerance) );
			}
			if( ++fragCS >= this.charge ) break;
		}
		
		Collections.sort(reporters);
		
		int start = 0;
		for( RawPeak rp : reporters ){
			for (int i=start; i<rawPL.size(); i++){
				if( rawPL.get(i).mz < rp.mz-rp.it ) continue;
				else if( rawPL.get(i).mz > rp.mz+rp.it ) {
					start = i;
					break;
				}			
				rawPL.remove(i);
				i--;
			}
		}	
	}
	
	public static ArrayList<MSMScan> getSpectra( String pList ) {
		
		String 	t_title="";
		int 	t_scanNo = 0, t_specIndex = 0;;
		double 	t_pmz= 0.;
		ArrayList<Integer> t_cslist = new ArrayList<Integer>();	

		StringTokenizer strSpec = new StringTokenizer(pList, "\r\n");
		if( strSpec.countTokens() < 3 || !strSpec.nextToken().startsWith("BEGIN IONS") )
			return new ArrayList<MSMScan>();
		
		String buf = "";
		while( strSpec.hasMoreTokens() ){
			buf = strSpec.nextToken();
			if( !buf.contains("=") ) break;
			
			if( buf.startsWith("TITLE=") ){
				t_title = buf.substring(buf.indexOf("=")+1).trim();	
				t_title = t_title.replaceAll("[/,:*?\"<>|\\\\]", "_");
			}
			else if( buf.startsWith("CHARGE=") ){												
				t_cslist = getCharge(buf.substring(buf.indexOf("=")+1));
			}			
			else if( buf.startsWith("PEPMASS=") ){						
				StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
				if( token.countTokens() != 0 ) t_pmz = Double.parseDouble( token.nextToken() );
			}
			else if( buf.startsWith("SCANS=") ){						
				StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
				if( token.countTokens() != 0 ) t_scanNo = Integer.parseInt( token.nextToken() );
			}
			else if( buf.startsWith("REINSEDCONDS=") ){						
			}
			else if( buf.startsWith("SPECINDEX=") ){						
				StringTokenizer token = new StringTokenizer(buf.substring(buf.indexOf("=")+1));
				if( token.countTokens() != 0 ) t_specIndex = Integer.parseInt( token.nextToken() );
			}
		}
		
		ArrayList<RawPeak> rawPL = new ArrayList<RawPeak>();				
		while( strSpec.hasMoreTokens() ){
			if( buf.startsWith("END") ) break;
			StringTokenizer token = new StringTokenizer(buf);
			if( token.countTokens() > 1 )
				rawPL.add( new RawPeak(Double.parseDouble(token.nextToken()), Double.parseDouble(token.nextToken())) );
			buf=strSpec.nextToken();
		}
		Collections.sort(rawPL);
		
		ArrayList<MSMScan> scanlist = new ArrayList<MSMScan>();
		MSMScan curScan = null;
		if( t_pmz != 0. ){
			if( t_cslist.size() > 0 ){
				for( Integer cs : t_cslist ){
					curScan = new MSMScan(t_title, t_specIndex, t_scanNo, t_pmz, cs);
					if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
				}
			}
			else {
				for(int cs=ScanIterator.MIN_ASSUMED_CHARGE; cs<=ScanIterator.MAX_ASSUMED_CHARGE; cs++){						
					curScan = new MSMScan(t_title, t_specIndex, t_scanNo, t_pmz, cs);
					if( curScan.setSpectrum(rawPL) ) scanlist.add(curScan);
				}
			}
		}

		return scanlist;
	}

	private static ArrayList<Integer> getCharge(String csline){
		ArrayList<Integer> cslist = new ArrayList<Integer>();
		
		StringTokenizer csTok = new StringTokenizer(csline, "|");
		while( csTok.hasMoreTokens() ){
			String csStr = csTok.nextToken();
			int t_cs = 0;
			int st = 0, ed = 1;
			for(int i=st; i<csStr.length(); i++){
				if( Character.isDigit( csStr.charAt(i) ) ){
					st = i;
					ed = st+1;
					break;
				}
			}
			for(int i=ed; i<csStr.length(); i++){
				if( !Character.isDigit( csStr.charAt(i) ) ){
					ed = i;
					break;
				}
			}
			try {
				t_cs= Integer.parseInt( csStr.substring(st,ed) );
			} catch (NumberFormatException e) {
				t_cs = 0;
			}
			if( t_cs!= 0 ) cslist.add(t_cs);
		}
		return cslist;
	}
	
}

























