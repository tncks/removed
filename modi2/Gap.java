package modi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import moda.ThreadLocalMutables;
import msutil.PGraph;

public class Gap implements SpecInterpretation {

	
	final Spectrum			sourceSpectrum;
	final Peptide				matchedPeptide;
	final int					start;
    final int end;		// both are inclusive
	final double				offset;
	final double				bStartMass;
	final double				yStartMass;		// since ver 1.0
	final PTMPosition			position;		// Nterm, C
	GapInterpretation	gapInterpretation;		
	// since ver 1.0
	boolean				interpreted;
	int					interpretID;
	
	final double				mostAbundantBYPeakIntensity;
	

	public				Gap(Peptide matchedPeptide, int start, int end, 
			double bStartMass, double yStartMass, double offset, 
			PTMPosition position, Spectrum sourceSpectrum, double mostAbundantBYPeakIntensity)
	{
		this.matchedPeptide = matchedPeptide;
		this.start = start;
		this.end = end;
		this.bStartMass = bStartMass;
		this.yStartMass = yStartMass;
		this.offset = offset;
		this.position = position;
		this.sourceSpectrum = sourceSpectrum;
		this.mostAbundantBYPeakIntensity = mostAbundantBYPeakIntensity;
	}
	public	void		setInterpretation( PTMSearchResult result )	
	{
		this.gapInterpretation = new GapInterpretation( this, result.getPTMRun() );
		this.interpretID = result.getID();
	}
	
	public	void		setInterpretation( PTMSearchResult result, PGraph graph )	
	{
		this.gapInterpretation = new GapInterpretation( this, result.getPTMRun(), graph );
		this.interpretID = result.getID();
	}
	
	public	void		setInterpreted(boolean interpreted) { this.interpreted = interpreted; }
	public	int			getInterpretID()	{ return interpretID; }
	public	boolean		isInterpreted()		{ return interpreted; }
	public	double		getOffset()			{ return offset; }
	public	double		getBStartMass()		{ return bStartMass; }
	public	double		getYStartMass()		{ return yStartMass; }
	public	Peptide		getMatchedPeptide() { return matchedPeptide; }
	public	Sequence	sequence() 			{ return matchedPeptide.subSequence(start, end+1); }
	public	int			getStart()			{ return start; }
	public	int			getEnd()			{ return end; }
	public	int			getLength()			{ return end - start + 1; }
	public	PTMPosition	getPosition()		{ return position; }
	public	GapInterpretation getGapInterpretation()	{ return gapInterpretation; }
	public	int			compareTo(SpecInterpretation t)	{ return this.start - t.getStart(); }
	


	public	int			getNumOfPTM()
	{
		if(gapInterpretation == null || gapInterpretation.isEmpty())
			return 0;
		else
			return gapInterpretation.getFirst().size();
	}
	@SuppressWarnings("StringBufferReplaceableByString")
    public	String		getSequenceStr()
	{ 
		StringBuffer output = new StringBuffer();
        output.append("_".repeat(Math.max(0, end - start + 1)));
		
		return output.toString();
	}
	// should be optimized later
	public	double		getScore()			
	{ 
		if(gapInterpretation == null || gapInterpretation.isEmpty())
			return 0.;
		else
			return gapInterpretation.getGapScore(0);
	}
	
	public	ArrayList<Peak> getBTheoreticalPeaks(PTMRun run)
	{
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> bPeaks = new ArrayList<>();
		double mass = bStartMass;
		for(int i=start; i<=end; i++)
		{
			mass = mass + matchedPeptide.get(i).getMass() + ptmMass[i-start];
			bPeaks.add(new Peak(-1, mass, 0.));
		}
		return bPeaks;
	}
	
	public	ArrayList<Peak> getYTheoreticalPeaks(PTMRun run)
	{
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> yPeaks = new ArrayList<>();
		double mass = bStartMass;
		for(int i=start; i<=end; i++)
		{
			yPeaks.add(new Peak(-1, sourceSpectrum.getCorrectedMW() - mass + 2*Constants.UNIT_MASS, 0.));
			mass = mass + matchedPeptide.get(i).getMass() + ptmMass[i-start];
		}
		return yPeaks;
	}
	
	public ArrayList<Peak> getTheoreticalPeaks(PTMRun run)
	{
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();

		ArrayList<Peak> allPeaks = new ArrayList<>();
		
		double mass = bStartMass, compMass= sourceSpectrum.getCorrectedMW() - bStartMass + 2*Constants.Proton;
		
		for(int i=start; i<end; i++) 
		{
			double addMass = matchedPeptide.get(i).getMass() + ptmMass[i-start];
			mass += addMass;
			compMass -= addMass;
			
			AminoAcid NTermResidue = matchedPeptide.get(i);
			AminoAcid CTermResidue = matchedPeptide.get(i+1);
			if(ptmMass[i-start] != 0)
				NTermResidue = AminoAcid.getAminoAcid(23);	// modified
			if(ptmMass[i-start+1] != 0)
				CTermResidue = AminoAcid.getAminoAcid(23);  // modified//*/
			

			{
				allPeaks.add(new MatchedPeak(-1, mass, 0, 1, PeakProperty.B_ION, NTermResidue, CTermResidue));
				if( sourceSpectrum.getCharge() > 2 ){
					allPeaks.add(new MatchedPeak(-1, (mass+1)/2, 0, 2, PeakProperty.B_ION, NTermResidue, CTermResidue));
				}
			}
			

			{
				allPeaks.add(new MatchedPeak(-1, compMass, 0, 1, PeakProperty.Y_ION, NTermResidue, CTermResidue));
				if( sourceSpectrum.getCharge() > 2 ){
					allPeaks.add(new MatchedPeak(-1, (compMass+1)/2, 0, 1, PeakProperty.Y_ION, NTermResidue, CTermResidue));
				}
			}
		}
		
		// add diagnostic peak
		HashSet<PTM> checkAdded = new HashSet<>();
		DiagnosticPTMIonDB diagPTMDB = new DiagnosticPTMIonDB();
		for(PTMOccurrence occr : run)
		{
			PTM ptm = occr.getPTM();
			if(!checkAdded.contains(ptm))
			{
				checkAdded.add(ptm);
				for(DiagnosticPTMIon ion : diagPTMDB.getDiagIon(ptm)){
					allPeaks.add( new DiagnosticPTMPeak(ion.getMass(), ion.getProb()) );
				}
			}
		}

		Collections.sort(allPeaks);
		
		return allPeaks;
	}

	public ArrayList<Peak> getTheoreticalPeaks()
	{
		if(gapInterpretation == null || gapInterpretation.isEmpty())
			return getTheoreticalPeaks(new PTMRun());
		else
			return getTheoreticalPeaks(gapInterpretation.getFirst());
	}
	
	public	static String getTheoreticalMassStr(ArrayList<Peak> peakList)
	{
		StringBuffer output = new StringBuffer();
		boolean start = true;
		for(Peak p : peakList)
		{
			if(start)
				start = false;
			else
				output.append(",");
			output.append(Constants.getString(p.getMass()));
		}
		return output.toString();
	}
	
	public double getScore(PTMRun run)
	{
		ArrayList<Peak> theoPeaks = this.getTheoreticalPeaks(run);
		
		double score = ( run.getError() > Mutables.fragmentTolerance )? -0.5 : 0;

		// inefficient
		int classification = FragmentationPattern.getClassification(sourceSpectrum, matchedPeptide);
		ArrayList<PeakPair> sharedPeaks = new PeakListComparator(sourceSpectrum, theoPeaks).getSharedPeaks();
		for(PeakPair p : sharedPeaks)
		{
			PeakProperty property = p.getSecond().getPeakProperty();
			double weight;
			if(property == PeakProperty.DIAGNOSTIC_PTM_ION)
				weight = ((DiagnosticPTMPeak)p.getSecond()).prob;
			else
				weight= 1;
			score += weight*p.getFirst().getNormIntensity();
		}
		
		score -= run.size()*0.5;//Constants.PTM_ADD_PENALTY;
	//	score -= run.getPenalty()*0.5;//Constants.PTM_ADD_PENALTY;
		
		return score;
	}

	boolean isInModifiedRange(double v) {
		if (Constants.minModifiedMass - (ThreadLocalMutables.get().gapTolerance) < v && v < Constants.maxModifiedMass + (ThreadLocalMutables.get().gapTolerance))
			return true;
		else return Math.abs(v) <= (ThreadLocalMutables.get().gapTolerance);
	}
	
	public boolean isReasonable(){
		if( !isInModifiedRange( this.offset ) ) return false;
		
	//	if ( this.offset > 0 ) return true;

        int pos = 0; //because no. of allowed mods is max two
		if( start == 0 ) pos = 1;
		if( end+1 == matchedPeptide.size() ) pos = 2;
		Sequence seq = matchedPeptide.subSequence(start, end+1);
		double minMod = Mutables.variableModifications.minimumModifiedMass(seq, pos);
		double maxMod = Mutables.variableModifications.maximumModifiedMass(seq, pos);
	//	System.out.println( this.offset + " " + maxMod + " " +minMod + " " + matchedPeptide.subSequence(start, end+1));
        return !(this.offset < (minMod - (ThreadLocalMutables.get().gapTolerance))) && !(this.offset > (maxMod + (ThreadLocalMutables.get().gapTolerance)));//*/
    }

	public double getScore(PTMRun run, PGraph graph)
	{
		double score = - run.getError();
		double [] ptmMass = new double[end-start+1];
		for(PTMOccurrence occr : run)
			ptmMass[occr.getPosition()] += occr.getPTM().getMassDifference();
		
		graph.refresh();
		ArrayList<Integer> singly = new ArrayList<>();
		ArrayList<Integer> doubly = new ArrayList<>();
		
		double mass = bStartMass, compMass= sourceSpectrum.getCorrectedMW() - bStartMass + 2*Constants.Proton;
		int count = 0;
		for(int i=start; i<end; i++)
		{
			double addMass = matchedPeptide.get(i).getMass() + ptmMass[i-start];
			mass += addMass;
			compMass -= addMass;
			count++;
			
			int index = graph.getIndexOfMatchedPeak( mass );
			if( index != -1 ){
				singly.add(index);
				score += graph.get(index).getNorm();
			}
			index = graph.getIndexOfMatchedPeak( compMass );
			if( index != -1 ){
				singly.add(index);
				score += graph.get(index).getNorm();
			}
			if( sourceSpectrum.getCharge() > 2 ){
				index = graph.getIndexOfMatchedPeak( (mass+Constants.Proton)/2 );
				if( index != -1 ){
					doubly.add(index);
					score += graph.get(index).getNorm();
				}
				index = graph.getIndexOfMatchedPeak( (compMass+Constants.Proton)/2 );
				if( index != -1 ){
					doubly.add(index);
					score += graph.get(index).getNorm();
				}
			}
		}
		for(int i=0; i<singly.size();i++){
			score += graph.getSupportingScore( singly.get(i++), 1, 1, 0.5 );
		}
		for(int i=0; i<doubly.size();i++){
			score += graph.getSupportingScore( doubly.get(i++), 1, 2, 0.5 );
		}
	//	score -= run.size()*2;
		score -= run.getPenalty()*2;
		return score;
	}
	
	
}