package modi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;

import moda.ThreadLocalMutables;
import modi.Constants;
import modi.Mutables;

import msutil.PGraph;
import msutil.PNode;

public class Spectrum extends ArrayList<Peak> implements Comparable<Spectrum> {

	private String		name = "";	// id
	private final double 		precursor;
	private final int 		charge;
	private final double 		observedMW;	// not used now
	private double 		correctedMW;
	private double		TIC= 0;
    private double		bpIntensity= 0;

    private	ArrayList<Peak> 	selectedPeak = null;
	
	public	Spectrum(double precursor, int charge) 
	{
		correctedMW = observedMW = (precursor - Constants.Proton)*charge;
		this.precursor = precursor;
		this.charge = charge;
	}
	public 	Spectrum(double precursor, int charge, String name)
	{
		this(precursor, charge);
		this.name = name;
	}
	public	void 	setExtraInformation(double bpi, double t){
		bpIntensity= bpi;
		TIC= t;
	}
	public	void 	addPeak(Peak p) 		{ this.add(p); }
	public	double	getObservedMW()	{ return observedMW; }
	public	double	getCorrectedMW()	{ return correctedMW; }
	public	void	setCorrectedParentMW(double mw)	{ 
		correctedMW = mw; 
	}
	
	public	double	getPrecursor()			{ return precursor; }
	public	double	getBasePaekIntensity()	{ return bpIntensity; }
	public	int		getCharge() 			{ return charge; }
	public	String	getName()				{ return name; }
	@SuppressWarnings("SameReturnValue")
    public	double	getRandomPeakIntensity() {
        return 0; }
	public  double  getTIC(){ return TIC; }
	@SuppressWarnings("SameReturnValue")
    public  double  getRTTime(){
        return 0; }
	public	Peak	getPeakByMass(double mass)
	{
		Peak result = new Peak(-1, -1, 0);
		for(Peak p : this)
		{
			if(Mutables.fEqual(mass, p.getMass()) && p.getIntensity() > result.getIntensity())
				result = p;
		}
		
		return result;
	}

	public	ArrayList<Peak>		getSelectedPeaks()	{ return selectedPeak; }
	
	public void normalizeIntensityLocally(){
		for(Peak p : this){			
			p.setNormIntensity(p.getIntensity()/getLocalMaxPeak(p.getMass()));	
		}
	}
	
	public double getLocalMaxPeak(double center){
		double lMax =0;
		int bin = 50;
		int index = binarySearchforPeaks( center-bin );
		
		while( this.get(index).getMass() <= center+bin ){
			if( this.get(index).getIntensity() > lMax )
				lMax= this.get(index).getIntensity();
			index++;
			if( index == this.size() )
				break;
		}	
		return lMax;
	}

	public PGraph getPeakGraph(){ //MOD SERIES

		int binSize= 100, considered= 10;
		double minimumCut= bpIntensity*0.001;
		double precursoeRange= (ThreadLocalMutables.get().precursorTolerance)/this.charge + Mutables.fragmentTolerance;
		
		PGraph graph = new PGraph(this.observedMW, charge);
		
		ArrayList<Peak> locBin = new ArrayList<>();
		
		int lowerBound = (int)this.get(0).getMass()+100, upperBound = (int)this.get(this.size()-1).getMass()-100;
		
		int startLocal = (int)this.get(0).getMass()/binSize;		
		int local= startLocal;
		
		for( Peak p : this ){
			if( local == (int)p.getMass()/binSize )
				locBin.add(p);
			else{
				locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
				
				int to= 1;
				for( Peak lp : locBin ){			
					if( startLocal == local && lp.intensity < minimumCut ) break;
					if( lp.mass < lowerBound || lp.mass > upperBound ) {
						if( lp.intensity < minimumCut ) {
							to++;
							if( considered < to ) break;
							continue;
						}//*/
					}
					
					if( this.precursor - lp.mass < precursoeRange && lp.mass - this.precursor < 2. ) {
						graph.add( new PNode(lp.mass, minimumCut, 10, Constants.rNorm[10]) );		
					}
					else //*/
					{
						graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );					
						to++;
					}
					if( considered < to ) break;
				}
				
				local = (int)p.getMass()/binSize;
				locBin.clear();
				locBin.add(p);
			}
		}

		locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
		int to= 1;
		for( Peak lp : locBin ){		
			if( lp.intensity < minimumCut ) break;				
			graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );
			to++;
			if( considered < to ) break;
		}//Last Region	

		graph.ready(TIC);
		
		return graph;
	}
	
	public PGraph getPeakGraph(int window, int pick, double noiseRatio){ // DBOND

        double minimumCut= bpIntensity*noiseRatio;
		
		PGraph graph = new PGraph( this.observedMW, charge );
		ArrayList<Peak> locBin = new ArrayList<>();

        int local= (int) this.get(0).getMass()/ window;
		for( Peak p : this ){
			if( local == (int)p.getMass()/ window)
				locBin.add(p);
			else{
				locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
				
				int to= 1;
				for( Peak lp : locBin ){			
					if( lp.intensity < minimumCut ) break;				
					graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );//Constants.rNorm[to]) );
					to++;
					if( pick < to ) break;
				}
				
				local = (int)p.getMass()/ window;
				locBin.clear();
				locBin.add(p);
			}
		}

		locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
		int to= 1;
		for( Peak lp : locBin ){		
			if( lp.intensity < minimumCut ) break;		
			graph.add( new PNode(lp.mass, lp.intensity, to, Constants.rNorm[to]) );
			to++;
			if( pick < to ) break;
		}//Last Region	
		graph.ready(TIC);
		return graph;
	}
	
	public PGraph getPeakGraph(int window, int pick) { // NOT used

        PGraph graph = new PGraph( this.observedMW, charge );
		ArrayList<Peak> locBin = new ArrayList<>();
					
		int local= (int)this.get(0).getMass()/ window;
		for( Peak p : this ){
			if( local == (int)p.getMass()/ window)
				locBin.add(p);
			else{
				locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
				
				int to= 1;
				for( Peak lp : locBin ){								
					graph.add( new PNode(lp.mass, lp.intensity, to, lp.normIntensity) );//Constants.rNorm[to]) );
					to++;
					if( pick < to ) break;
				}				
				local = (int)p.getMass()/ window;
				locBin.clear();
				locBin.add(p);
			}
		}

		locBin.sort(Collections.reverseOrder( new IntensityComparator()) );
		int to= 1;
		for( Peak lp : locBin ){							
			graph.add( new PNode(lp.mass, lp.intensity, to, lp.normIntensity) );
			to++;
			if( pick < to ) break;
		}//Last Region	
		graph.ready(TIC);
		return graph;
	}
	
	public 	String	toString() 
	{
		StringBuffer buffer = new StringBuffer();
		buffer.append("Spectrum("+precursor+","+charge+")\n");
		Iterator<Peak> e = this.iterator();
		while(e.hasNext())
		{
			buffer.append(e.next().toString());
		}
		return buffer.toString();
	}
	
	public int peakSelection (double selectionWindowSize, int numOfPeaksInWindow ){
		int globalSelectedSize = (int) (this.observedMW / selectionWindowSize * numOfPeaksInWindow);
		return peakSelection(globalSelectedSize, selectionWindowSize, numOfPeaksInWindow);
	}
	public int peakSelection (int globalSelectedSize, double selectionWindowSize, int numOfPeaksInWindow ){
		assert(this.checkSorted());
		if(globalSelectedSize > this.size())
			globalSelectedSize = this.size();
		if(this.size() == 0)
			return -1;	// later : error code
		// global selection
		Iterator<Peak> it = this.iterator();
		TreeSet<Peak> globalSelected = new TreeSet<>(new IntensityComparator());
		for(int i=0; it.hasNext() && i<globalSelectedSize; i++)
			globalSelected.add(it.next());
		TreeSet<Peak> notSelected = new TreeSet<>(new MassComparator());
		Peak curPeak;
		while(it.hasNext()) 
		{
			curPeak = it.next();
			if(curPeak.getIntensity() > globalSelected.first().getIntensity())
			{
				notSelected.add(globalSelected.first());
				globalSelected.remove(globalSelected.first());	// delete smallest peak
				globalSelected.add(curPeak);
			}
			else
				notSelected.add(curPeak);
		}
		
		// local selection
		TreeSet<Peak> selected = new TreeSet<>(new MassComparator());
		selected.addAll(globalSelected);
		int currentBinPeakSize;
		double currentBinStart, currentBinEnd;
		// do not select additional peaks in (0, selectionWindowSize)
		for(int i=2; i<(int)(correctedMW/selectionWindowSize+1)*2; i++)
		{
			currentBinStart = i*selectionWindowSize/2;
			currentBinEnd = currentBinStart + selectionWindowSize;
			currentBinPeakSize = selected.subSet(new Peak(0, currentBinStart, 0), new Peak(0, currentBinEnd, 0)).size();
			if( currentBinPeakSize < numOfPeaksInWindow ){
				Peak [] newPeaks = notSelected.subSet(new Peak(0, currentBinStart, 0.), new Peak(0, currentBinEnd, 0.)).toArray(new Peak[0]);
				Arrays.sort(newPeaks, Collections.reverseOrder(new IntensityComparator()));
				for(int j=0; j<newPeaks.length && j < numOfPeaksInWindow - currentBinPeakSize; j++)
				{
					if( newPeaks[j].getNormIntensity() > Constants.minNormIntensity ){
						notSelected.remove(newPeaks[j]);
						selected.add(newPeaks[j]);
					}
					else break;
				} 
			}
		}
		
		selected.add(new Peak(-1, Constants.B_ION_OFFSET+Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_B_ION_ONLY));
		selected.add(new Peak(-1, Constants.Y_ION_OFFSET+Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_Y_ION_ONLY));
		selected.add(new Peak(-1, correctedMW-Constants.H2O+Constants.Proton-Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_B_ION_ONLY));
		selected.add(new Peak(-1, correctedMW+Constants.Proton-Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_Y_ION_ONLY));

		selectedPeak = new ArrayList<>(selected);
		setScoreOfSelectedPeaks(selectedPeak, 1, Mutables.massToleranceForDenovo);

        return selectedPeak.size();
	}	

	public 	void	printSelectedPeak()
	{
		System.out.println(selectedPeak.size() + " peaks are selected");
		for(Peak peak : selectedPeak)
			System.out.println(peak);
	}

	double getMassDifference(Peak p1, Peak p2)
	{
		return Math.abs(p1.mass - p2.mass);
	}

	boolean		extendable(Tag t1, Tag t2)
	{
		boolean ex= t1.getLast() == t2.getFirst();
		return ex;
	}

	TagPool extendTags(TagPool target, TagPool source)
	{
		TagPool result = new TagPool();
		for(int i=0; i<target.size(); i++)
			for(int j=0; j<source.size(); j++)
				if(extendable(target.get(i), source.get(j)))
					result.add(Tag.merge(target.get(i), source.get(j)));

		return result;
	}
	
	public	TagPool	generateTags( int minTagLength, int minTagLengthPeptideShouldContain, double massTolerance ) {// by NA
		// selectedPeak array should be sorted before generating Tags
		assert(selectedPeak != null);
	//	double maxGap = AminoAcid.getAminoAcid('W').getMonoMass() + massTolerance; // Mass of 'W' is biggest
		double maxGap = AminoAcid.getMaxMonoMass() + massTolerance; // Mass of 'W' is biggest
		
		TagPool tags = new TagPool();		
		for(int i=0; i<selectedPeak.size()-1; i++)
		{		
			for(int j=i+1; j<selectedPeak.size(); j++){
				double diff = getMassDifference(selectedPeak.get(i), selectedPeak.get(j));
				if( diff > maxGap )
					break;
				
				// one length tag
				ArrayList<AminoAcid> codeList = AminoAcid.getCode( diff, massTolerance );//massTolerance);
				for(AminoAcid code : codeList) {
					tags.add(new Tag(selectedPeak.get(i), selectedPeak.get(j), code, this));
				}							
			}			
		}
		
		TagPool oneLengthTagPool 	= (TagPool)tags.clone();
		TagPool tempTags			= (TagPool)tags.clone();
		
		for(int i=1; i<minTagLengthPeptideShouldContain; i++)
		{
			tempTags = extendTags(tempTags, oneLengthTagPool);
			tags.addAll(tempTags);
		}//*/		
		
		if( this.charge > 2 ){
			tags.addAll( this.generateDoublyTags(minTagLength, minTagLengthPeptideShouldContain, massTolerance*2) );			
		}//*/
		
		tags.setTagScores();
		tags = tags.extractQualifiedTag(minTagLength, Constants.MAX_TAG_SIZE);//*/	

        return tags;
	}
	
	private	TagPool	generateDoublyTags(int minTagLength, int minTagLengthPeptideShouldContain, double massTolerance)
	{// by NA
		// selectedPeak array should be sorted before generating Tags
		assert(selectedPeak != null);
		
		ArrayList<Peak> virtualDoublyPeak = new ArrayList<>();
		for(int i=0; i<selectedPeak.size(); i++){
			
			if( selectedPeak.get(i).property == PeakProperty.N_TERM_B_ION_ONLY || 
					selectedPeak.get(i).property == PeakProperty.N_TERM_Y_ION_ONLY ){
				virtualDoublyPeak.add( selectedPeak.get(i) );
				continue;
			}
			
			if( selectedPeak.get(i).property == PeakProperty.C_TERM_B_ION_ONLY || 
					selectedPeak.get(i).property == PeakProperty.C_TERM_Y_ION_ONLY ){
				virtualDoublyPeak.add( selectedPeak.get(i) );
				continue;
			}
			
			if( selectedPeak.get(i).getMass()*2 < correctedMW )
			{
				Peak virtualPeak = new Peak(-1, selectedPeak.get(i).getMass()*2-Constants.Proton, selectedPeak.get(i).getIntensity(), 1, PeakProperty.VIRTUAL_PEAK);
				virtualPeak.setNormIntensity( selectedPeak.get(i).getNormIntensity() );				
				virtualDoublyPeak.add( virtualPeak );
			}
		}
		Collections.sort( virtualDoublyPeak );
		setScoreOfSelectedPeaks(virtualDoublyPeak, 2, massTolerance);
		
	//	double maxGap = AminoAcid.getAminoAcid('W').getMass() + massTolerance; // Mass of 'W' is biggest
		double maxGap = AminoAcid.getMaxMonoMass() + massTolerance; // Mass of 'W' is biggest
		TagPool tags = new TagPool();		
		for(int i=0; i<virtualDoublyPeak.size()-1; i++)
		{
			for(int j=i+1; j<virtualDoublyPeak.size(); j++){
				double diff = getMassDifference(virtualDoublyPeak.get(i), virtualDoublyPeak.get(j));
				if( diff > maxGap )
					break;
				
				// one length tag
				ArrayList<AminoAcid> codeList = AminoAcid.getCode(diff, massTolerance);
				for(AminoAcid code : codeList){//public Peak(int index, double mass, double intensity, int charge, PeakProperty property)				
					tags.add(new Tag(virtualDoublyPeak.get(i), virtualDoublyPeak.get(j), code, this));
				}							
			}			
		}
		
		TagPool oneLengthTagPool 	= (TagPool)tags.clone();
		TagPool tempTags			= (TagPool)tags.clone();
		
		for(int i=1; i<minTagLengthPeptideShouldContain; i++) {
			tempTags = extendTags(tempTags, oneLengthTagPool);
			tags.addAll(tempTags);
		}	

		return tags;
	}
	
	public Tag getB2Tag(Sequence seq)
	{
		assert(selectedPeak != null);
		
		int virB0 = -1;
		for(int j=0; j<selectedPeak.size(); j++){
			if( selectedPeak.get(j).getPeakProperty() == PeakProperty.N_TERM_B_ION_ONLY ){
				virB0 = j;
				break;
			}
		}
		
		double b2mz = selectedPeak.get(virB0).getMass() + seq.getMonoMass();
		
		for(int j=virB0; j<selectedPeak.size(); j++){
			
			double diff = b2mz - selectedPeak.get(j).getMass();			
			if( diff < -Mutables.fragmentTolerance ) break;
			
			if( Math.abs(diff) <= Mutables.fragmentTolerance ){
			//	System.out.println(name + " | " + b2mz);
				Peak virtualPeak = new Peak(-1, b2mz - seq.getMonoMass(1, 2), 0, 1, PeakProperty.VIRTUAL_PEAK);
				virtualPeak.setNormIntensity(0);
				virtualPeak.setProbability(0);
				Tag t1 = new Tag(selectedPeak.get(virB0), virtualPeak, seq.get(0), this);
				Tag t2 = new Tag(virtualPeak, selectedPeak.get(j), seq.get(1), this);				
				return Tag.merge(t1, t2);	
			}
		}//*/
		return null;
	}
	
	public int binarySearchforPeaks( double left )
	{
		int index;	
		if( left <= this.get(0).getMass() )
			index= 0;
		else if( left > this.get(this.size()-1).getMass() )
			index= this.size()-1;
		else
		{
			int M, L= 0, R= this.size()-1;
			while( R - L > 1 )
			{
				M= ( L + R ) /2;

				if( left <= this.get(M).getMass() )
					R= M;
				else
					L= M;
			}
			index= R;
		}	
		return index;
	}
	
	public double getMatchedPeak( double mz ){
		double it=0;

		int id= binarySearchforPeaks( mz- Mutables.fragmentTolerance );
		if( this.get(id).getMass() < mz - Mutables.fragmentTolerance )
			return it;
		
		while( this.get(id).getMass() <= mz + Mutables.fragmentTolerance )
		{
			if( this.get(id).getNormIntensity() > it )
				it = this.get(id).getNormIntensity();		
			id++;
			if( id == this.size() )
				break;
		}
		return it;
	}
	
	public double getPeakEvidence( Peak p ){
		double mz= p.getMass();	
		double ev= p.getNormIntensity(), it;
		
		if( p.property == PeakProperty.N_TERM_B_ION_ONLY || 
				p.property == PeakProperty.N_TERM_Y_ION_ONLY)
			return 1;
		
		if( p.property == PeakProperty.C_TERM_B_ION_ONLY || 
				p.property == PeakProperty.C_TERM_Y_ION_ONLY)
			return 1;
	
		if( p.getNormIntensity() < (it = getMatchedPeak(mz-Constants.IsotopeSpace)) ){ return ev; }	//Isotope	
		
		ev += getMatchedPeak( mz+1.0 );
		if( p.getNormIntensity() >= (it=getMatchedPeak(mz-Constants.NH3)) ){ ev += it; }	//NH3
		if( p.getNormIntensity() >= (it=getMatchedPeak(mz-Constants.H2O)) ){ ev += it; }	//H2O

		ev += getMatchedPeak( p.getComplementMass(correctedMW) ); 		//Complementary */
		
		return ev;
	}
	
	public boolean isConfidentPeak( Peak tar, double factor ){
	//	System.out.println(tar.mass +" "+ tar.getIndex() + " " + this.size());
		int me = tar.getIndex();
		if( me < 0 ) return true;
		if( me > 0  && tar.getMass() - this.get(me-1).getMass() < 1.5 ){
			if( tar.getIntensity()/5 < this.get(me-1).getIntensity() ) return false;
		}

		if( me < this.size()-1 && this.get(me+1).getMass() - tar.getMass() < 1.5 ){
            return !(tar.getIntensity() < this.get(me + 1).getIntensity());
		}			
			
		return true;		
	}

	
	public boolean checkSorted()	// identical to Tag's method
	{
		Peak tmp = null;
		for (Peak p : this)
		{
			if (tmp != null)
			{
				if (p.compareTo(tmp)<0) return false;
			}
			tmp = p;
		}
		return true;
	}
	
	public int compareTo(Spectrum s) 
	{
		if( observedMW > s.observedMW ) return 1;
		else if( observedMW == s.observedMW ) return 0;
		else return -1;
	}
	
	public String toDTA(){
		
		StringBuffer buffer = new StringBuffer();
		buffer.append( (this.observedMW+Constants.Proton) + "\t" + this.getCharge()+"\r\n");
		
		for(Peak p : this){
			buffer.append(p.getMass() + "\t" + p.getIntensity()+"\r\n");
		}
		buffer.append("\r\n");			
		return buffer.toString();
}
	
	private void setScoreOfSelectedPeaks( ArrayList<Peak> selected, int assumedCS, double tolerance )
	{
		double isotopeDelta = Constants.IsotopeSpace, NH3Delta = Constants.NH3, H2ODelta = Constants.H2O;
		for( int node=0; node<selected.size(); node++){
			
			if( selected.get(node).property == PeakProperty.N_TERM_B_ION_ONLY || 
					selected.get(node).property == PeakProperty.N_TERM_Y_ION_ONLY ){
				selected.get(node).setProbability(1.);
				continue;
			}
			
			if( selected.get(node).property == PeakProperty.C_TERM_B_ION_ONLY || 
					selected.get(node).property == PeakProperty.C_TERM_Y_ION_ONLY ){
				selected.get(node).setProbability(1.);
				continue;
			}
			
			double targetmz;
			double score = selected.get(node).getNormIntensity(); 
			double ionMZ= selected.get(node).getMass(); 
			double ionIT= selected.get(node).getIntensity();
			
			// check this is isotope???
			int OK = -1;    	
			double okmax = 0;
			targetmz= ionMZ-isotopeDelta;
			for(int i=node-1; i>-1; i--){
				if( selected.get(i).getMass() < targetmz-tolerance ) break;			
				else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
					if( selected.get(i).getIntensity() > okmax ) {
						okmax = selected.get(i).getIntensity();
						OK = i;
					}
				}
			}
			if( OK != -1 ) {
				if( okmax > ionIT ) selected.get(node).setProbability(0.);
				else selected.get(node).setProbability(score);
				continue;
			}//*/
			

			int ISO = node + 1; // plus isotope peak    
			double prevISO = ionIT;
			boolean isoDecent = false;
			targetmz= ionMZ + isotopeDelta;
			for( int nth_iso=1 ;  ; nth_iso++  ){
				int H = -1;
				double imax = 0;				
				for(int i=ISO; i<selected.size(); i++){
					if( selected.get(i).getMass() > targetmz+Mutables.massToleranceForDenovo ) break;
					else if( Math.abs(selected.get(i).getMass()-targetmz) < Mutables.massToleranceForDenovo ){
						if( selected.get(i).getIntensity() > imax ) {
							imax = selected.get(i).getIntensity();
							H = i;
						}
					}
				}

				if( H == -1 || prevISO < imax/5 ) break;
				if( isoDecent && prevISO < imax ) break;
				if( prevISO > imax ) isoDecent=true;
				
				score += selected.get(H).getNormIntensity()/nth_iso;	
				ISO = H+1;
				targetmz= selected.get(H).getMass() + isotopeDelta;		
				prevISO = imax;
			}//*/

			int NLOSS = -1;    	
			double lossmax=0, lossmz=0;
			targetmz= ionMZ-H2ODelta;
			for(int i=node-1; i>-1; i--){
				if( selected.get(i).getMass() < targetmz-tolerance ) break;			
				else if( Math.abs(selected.get(i).getMass()-(ionMZ-NH3Delta)) < tolerance ){
					if( selected.get(i).getIntensity() > lossmax ) {
						lossmax = selected.get(i).getIntensity();
						lossmz = selected.get(i).getMass();
						NLOSS = i;
					}
				}
				else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
					if( selected.get(i).getIntensity() > lossmax ) {
						lossmax = selected.get(i).getIntensity();
						lossmz = selected.get(i).getMass();
						NLOSS = i;
					}
				}
			}
			if( NLOSS != -1 ){
				double lossScore = selected.get(NLOSS).getNormIntensity();				
				int NL_H = -1; // plus isotope peak    	
				double nsmax = 0;
				targetmz= lossmz+isotopeDelta;
				for(int i=NLOSS+1; i<selected.size(); i++){
					if( selected.get(i).getMass() > targetmz+tolerance ) break;			
					else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
						if( selected.get(i).getIntensity() > nsmax ) {
							nsmax = selected.get(i).getIntensity();
							NL_H = i;
						}
					}
				}
				if( NL_H != -1 ) lossScore += selected.get(NL_H).getNormIntensity();				
				score += lossScore*0.5;
				
			}//*/
		
			double complement = correctedMW - ionMZ + Constants.Proton*2;	
			score += getMatchedPeak( complement ); 
			if( assumedCS == 1 && charge > 2 ) {
				score += getMatchedPeak( (complement+Constants.Proton)/2 ); 
			}			
			selected.get(node).setProbability(score);
		}
	}
}

class PrecursorComparator implements Comparator<Spectrum>
{
	public int compare(Spectrum spec1, Spectrum spec2)
	{
		double preMass1 = spec1.getPrecursor();
		double preMass2 = spec2.getPrecursor();
		if(preMass1 > preMass2)
			return 1;
		else if(preMass1 == preMass2)
			return 0;
		else 
			return -1;
	}
	
	public boolean equals(Spectrum spec1, Spectrum spec2)
	{
		return spec1.getPrecursor() == spec2.getPrecursor();
	}
}

class SpecNameComparator implements Comparator<Spectrum>
{
	public int compare(Spectrum spec1, Spectrum spec2)
	{	
		if( spec1.getName().compareTo(spec2.getName()) > 0 )
			return 1;
		else if( spec1.getName().compareTo(spec2.getName()) == 0 )
			return 0;
		else 
			return -1;
	}
	
	public boolean equals(Spectrum spec1, Spectrum spec2)
	{
		return spec1.getName().compareTo(spec2.getName()) == 0;
	}
}
