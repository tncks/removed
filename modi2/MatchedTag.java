package modi;

import java.util.ArrayList;
import java.util.Collections;

public class MatchedTag extends Tag implements SpecInterpretation {
	Peptide			matchedPeptide;
	int				start;
	int				end;
	IonDirection	direction;
	double			score = 0;


	@SuppressWarnings("CopyConstructorMissesField")
    public	MatchedTag(MatchedTag tag)
	{
		super(tag);
		this.matchedPeptide = tag.matchedPeptide;
		this.start = tag.getStart();
		this.end = tag.getEnd();
		this.direction = tag.direction;
	}
	
	public	MatchedTag(Tag tag, PeptideDBHit hitInformation, IonDirection direction)
	{
		super(tag);
		this.matchedPeptide = hitInformation.getHitPeptide();
		this.start = hitInformation.getStart();
		this.end = hitInformation.getEnd();
		tagSequence = matchedPeptide.subSequence(start, end+1);
		this.direction = direction;		
	}

	public	Peptide			getMatchedPeptide()	{ return matchedPeptide; }
	public	int				getStart()			{ return start; }
	public	int				getEnd()			{ return end; }
	public	IonDirection	getDirection() 		{ return direction; }
	public	int				compareTo(SpecInterpretation t)	{ return this.start - t.getStart(); }	
	public	double			getNTermOffset()	
	{ 
		if( direction == IonDirection.B_DIRECTION )
			return getBIonNtermOffset() - matchedPeptide.getMonoMass(0, start);
		else
			return getYIonNtermOffset() - matchedPeptide.getMonoMass(0, start);
	}
	public	double			getNTermOffsetByLastPeak()	
	{ 
		if(direction == IonDirection.B_DIRECTION)
			return this.get(this.size()-1).getMass() - Constants.B_ION_OFFSET - matchedPeptide.getMonoMass(0, end+1);
		else
			return this.getSourceSpectrum().getCorrectedMW() - this.getFirst().getMass() + Constants.UNIT_MASS - matchedPeptide.getMonoMass(0, end+1);
	}

	public	double			getCTermOffset()	
	{ 
		if( direction == IonDirection.B_DIRECTION )
			return getBIonCtermOffset() - matchedPeptide.getMonoMass(end+1, matchedPeptide.size());
		else
			return getYIonCtermOffset() - matchedPeptide.getMonoMass(end+1, matchedPeptide.size());
	}
	public	double			getOffset()			
	{ 
		if(direction == IonDirection.B_DIRECTION)
			return getNTermOffset();
		else
			return getCTermOffset();
	}
	private AminoAcid		getNTermResidue(int peakIndex)
	{
		if(this.direction == IonDirection.B_DIRECTION)
			return this.matchedPeptide.get(start+peakIndex-1);
		else
			return this.matchedPeptide.get(start+peakIndex);
	}
	private AminoAcid		getCTermResidue(int peakIndex)
	{
		if(peakIndex < 0 && peakIndex >= this.matchedPeptide.size())
			return null;
		if(this.direction == IonDirection.B_DIRECTION)
			return this.matchedPeptide.get(start+peakIndex);
		else
			return this.matchedPeptide.get(start+peakIndex-1);
	}
	public	void			setScore()
	{		
		score= 0;
		for(Peak p: this){
			score += p.getProbability();
		}					
	}
	
	public	double			getScore()			{ return score; }
	public	String			getSequenceStr()	{ return this.tagSequence.toString(direction); } 

	public boolean checkCompatibility()
	{
		if( start == 0 && !Mutables.fEqual(getNTermOffset()-Constants.NTERM_FIX_MOD, 0) )
			return false;
        return end != matchedPeptide.size() - 1 || Mutables.fEqual(getCTermOffset() - Constants.CTERM_FIX_MOD, 0);
    }
	


	//	 called when this and target are overlapped or adjacent while matchedPeptide, direction, offsets are same
	public	void	extend(MatchedTag target)	 
	{
		assert(this.matchedPeptide==target.matchedPeptide
				&& this.direction == target.direction
				&& Mutables.fEqual( this.getOffset(), target.getOffset()) );
		
		RelativePosition ir= this.getRelativePosition(target);
		if( ir == RelativePosition.SEPERATED )
			return;

		start = Math.min(getStart(), target.getStart());
		end = Math.max(getEnd(), target.getEnd());

		// Peak merge : should be optimized later
		PeakListComparator.mergePeakList(this, target);
		
		this.tagSequence = matchedPeptide.subSequence(start, end+1);
	}

	public RelativePosition getRelativePosition(MatchedTag tag)
	{
		if(this.getStart() <= tag.getStart() && this.getEnd() >= tag.getEnd())
			return RelativePosition.INCLUDING;
		else if(tag.getStart() <= this.getStart() && tag.getEnd() >= this.getEnd())
			return RelativePosition.INCLUDED;
		else if((tag.getStart() > this.getStart() && tag.getStart() <= this.getEnd()) || // i and j overlap
			(this.getStart() > tag.getStart() && this.getStart() <= tag.getEnd()))
			return RelativePosition.OVERLAP;	// this and tag are overlaped
		else if((this.getEnd()+1 == tag.getStart()) || // i and j adjacent
				(tag.getEnd()+1 == this.getStart()))
			return RelativePosition.ADJACENT;	// this and tag are adjacent
		else
			return RelativePosition.SEPERATED;
	}
	
	public ArrayList<Peak> getTheoreticalPeaks()
	{
		ArrayList<Peak> allPeaks = new ArrayList<>();

		if(this.direction == IonDirection.B_DIRECTION)
		{
		
			for(int i=0; i<this.size(); i++)
			{
				Peak p = this.get(i);
				double mass = p.getMass();
				double compMass = p.getComplementMass(sourceSpectrum.getCorrectedMW());
				allPeaks.add(new MatchedPeak(-1, mass, 0, 1, PeakProperty.B_ION, getNTermResidue(i), getCTermResidue(i)));
				allPeaks.add(new Peak(-1, mass-17*Constants.UNIT_MASS, 0, 1, PeakProperty.B_MINUS_NH3_ION));
				allPeaks.add(new Peak(-1, mass+Constants.UNIT_MASS, 0, 1, PeakProperty.ISOTOPE));
				allPeaks.add(new Peak(-1, mass-18*Constants.UNIT_MASS, 0, 1, PeakProperty.B_MINUS_H2O_ION));
				allPeaks.add(new Peak(-1, mass-28*Constants.UNIT_MASS, 0, 1, PeakProperty.A_ION));

				allPeaks.add(new MatchedPeak(-1, compMass, 0, 1, PeakProperty.Y_ION, getNTermResidue(i), getCTermResidue(i)));
				allPeaks.add(new Peak(-1, compMass+Constants.UNIT_MASS, 0, 1, PeakProperty.ISOTOPE));
				allPeaks.add(new Peak(-1, compMass-17*Constants.UNIT_MASS, 0, 1, PeakProperty.Y_MINUS_NH3_ION));
				allPeaks.add(new Peak(-1, compMass-18*Constants.UNIT_MASS, 0, 1, PeakProperty.Y_MINUS_H2O_ION));
			}
		}
		else
		{
			for(int i=0; i<this.size(); i++)
			{
				Peak p = this.get(i);
				double mass = p.getMass();
				double compMass = p.getComplementMass(sourceSpectrum.getCorrectedMW());

				allPeaks.add(new MatchedPeak(-1, compMass, 0, 1, PeakProperty.B_ION, getNTermResidue(i), getCTermResidue(i)));
				allPeaks.add(new Peak(-1, compMass+Constants.UNIT_MASS, 0, 1, PeakProperty.ISOTOPE));
				allPeaks.add(new Peak(-1, compMass-17*Constants.UNIT_MASS, 0, 1, PeakProperty.B_MINUS_NH3_ION));
				allPeaks.add(new Peak(-1, compMass-18*Constants.UNIT_MASS, 0, 1, PeakProperty.B_MINUS_H2O_ION));
				allPeaks.add(new Peak(-1, compMass-28*Constants.UNIT_MASS, 0, 1, PeakProperty.A_ION));

				allPeaks.add(new MatchedPeak(-1, mass, 0, 1, PeakProperty.Y_ION, getNTermResidue(i), getCTermResidue(i)));
				allPeaks.add(new Peak(-1, mass+Constants.UNIT_MASS, 0, 1, PeakProperty.ISOTOPE));
				allPeaks.add(new Peak(-1, mass-17*Constants.UNIT_MASS, 0, 1, PeakProperty.Y_MINUS_NH3_ION));
				allPeaks.add(new Peak(-1, mass-18*Constants.UNIT_MASS, 0, 1, PeakProperty.Y_MINUS_H2O_ION));
			}
		}
		Collections.sort(allPeaks);
		return allPeaks;
	}

	public	String toString()
	{
		return tagSequence.toString(direction) + " " + super.toString();
	}	
	
	public	String getSequence()
	{
		return tagSequence.toString(direction);
	}

	public void trim(int start, int len) {
		int i;
		if( this.direction == IonDirection.B_DIRECTION ){
			if( this.start==start ){
				for(i=0;i<len;i++){
					this.removeFirst();
					this.tagSequence.removeFirst();
					this.start++;				
				}
			}
			else{
				for(i=0;i<len;i++){
					this.remove(this.size()-1);
					this.tagSequence.removeLast();
					this.end--;
				}
			}
		}
		else{
			if( this.start==start ){
				for(i=0;i<len;i++){
					this.remove(this.size()-1);
					this.start++;
					this.tagSequence.removeFirst();
				}			
			}
			else{
				for(i=0;i<len;i++){
					this.removeFirst();
					this.tagSequence.removeLast();
					this.end--;
				}
			}
		}
	}
	public MatchedTag getTrimedTag(int start, int len) {
		if( this.sequence().size() == len )
			return null;
		int i;
		MatchedTag tag= new MatchedTag(this);
		if( tag.direction == IonDirection.B_DIRECTION ){
			if( tag.start==start ){
				for(i=0;i<len;i++){
					tag.removeFirst();
					tag.tagSequence.removeFirst();
					tag.start++;				
				}
			}
			else{
				for(i=0;i<len;i++){
					tag.removeLast();
					tag.tagSequence.removeLast();
					tag.end--;
				}
			}
		}
		else{
			if( this.start==start ){
				for(i=0;i<len;i++){
					tag.removeLast();
					tag.start++;
					tag.tagSequence.removeFirst();
				}			
			}
			else{
				for(i=0;i<len;i++){
					tag.removeFirst();
					tag.tagSequence.removeLast();
					tag.end--;
				}
			}
		}
		return tag;
	}
	
	public boolean isNormalSame(MatchedTag tag){
		if( this.direction != tag.direction )
			return false;
		if( !this.tagSequence.toString().equals(tag.tagSequence.toString()) )
			return false;
        return Mutables.fEqual(this.getOffset(), tag.getOffset());
    }
	
	public boolean isComplementarySame(MatchedTag tag){
		if( this.direction == tag.direction )
			return false;
		if( !this.tagSequence.toString().equals(tag.tagSequence.toString()) )
			return false;
        return Mutables.fEqual(this.getNTermOffset(), tag.getNTermOffset());
    }
	
	public boolean isLikelyChild(MatchedTag tag){
		if( this.direction != tag.direction )
			return false;

		if( Mutables.fEqual( this.getNTermOffset()-tag.getNTermOffset(), Constants.IsotopeSpace) )
			return true;
        return Mutables.fEqual(this.getNTermOffset() - tag.getNTermOffset(), -Constants.IsotopeSpace);
    }

	public boolean useSamePeak(MatchedTag tag){
		
	//	if( this.direction != tag.direction && this.sequence().toString().compareTo(tag.sequence().toString())!=0 ) return false;	
		if( this.direction != tag.direction ) return false;	
		
		for( Peak p : this ){
			for( Peak t : tag ){
				if( p.equals(t) ) return true;
			}
		}
		
		return false;
	}
	
	@SuppressWarnings("RedundantMethodOverride")
    public int hashCode(){
		return tagSequence.toString().hashCode();
	}
	public boolean equals(Object o){		
		if(!(o instanceof MatchedTag tag))
			   return false;

        return ( isComplementarySame(tag) || isNormalSame(tag) );
	}
	public	MatchedTag(Tag tag)
	{
		super(tag);
	}
	public static MatchedTag getMatchedB2Tag(Tag tag, Peptide pept)
	{
	    MatchedTag b2m = new MatchedTag(tag);
	    b2m.matchedPeptide = pept;
	    b2m.start = 0;
	    b2m.end = 1;
	    b2m.tagSequence = pept.subSequence(0, 2);
		b2m.direction = IonDirection.B_DIRECTION;
		return b2m;
	}
}
