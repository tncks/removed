package modi;

import java.util.ArrayList;

public class Peptide extends Sequence {
		
	private ArrayList<SrcProteinInfo> srcProteinInfo = new ArrayList<>();
	
	public Peptide( int srcProteinID, int start, int end ) {
		srcProteinInfo.add( new SrcProteinInfo( srcProteinID, start, end ));
	}
	
	public Peptide( String sequence, ArrayList<SrcProteinInfo> src) {
		this.addAll( Sequence.getSequence(sequence) );
		this.srcProteinInfo = src;
	}
	
	public Peptide( Peptide peptide ) {
		this.addAll( peptide );
		this.srcProteinInfo = peptide.srcProteinInfo;
	}
	
	public ArrayList<SrcProteinInfo> getSrcProteinInfo() {
		return srcProteinInfo;
	}
	public int getNTT(){
		return srcProteinInfo.getFirst().getNTT();
	}
	public int getStartCoordinate(){
		return srcProteinInfo.getFirst().getStartPos();
	}
	public int getEndCoordinate(){
		return srcProteinInfo.getFirst().getEndPos();
	}

	public double getMolecularMass()
	{
		return this.getMonoMass() + Constants.H2O;
	}

	public	int compareTo( Peptide p ) {
		for (int i=0; i<size() && i<p.size(); i++)
		{
			int r = get(i).compareTo(p.get(i));
			if (r<0) return -1;
			else if (r>0) return 1;
		}
		int r = size()-p.size();
		if (r<0) return -1;
		else if (r>0) return 1;
		else return 0;		
	}

	public int hashCode(){
		return srcProteinInfo.hashCode();
	}
	
	public boolean equals(Object o){		
		if(!(o instanceof Peptide p))
			   return false;
        return this.srcProteinInfo.hashCode() == p.srcProteinInfo.hashCode();
	}
	
}
