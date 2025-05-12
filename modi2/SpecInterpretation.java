package modi;

import java.util.ArrayList;
import java.util.Comparator;

public interface SpecInterpretation extends Comparable<SpecInterpretation> {
	Sequence	sequence();
	Peptide		getMatchedPeptide();
	int			getStart();
	int			getEnd();
	String		getSequenceStr();
	double		getScore();
	int			compareTo(SpecInterpretation t);
	ArrayList<Peak> getTheoreticalPeaks();

}

class SpecInterpretationComparator implements Comparator<SpecInterpretation> {
	public	int	compare(SpecInterpretation t1, SpecInterpretation t2)
	{
		return t1.getStart() - t2.getStart();
	}
	
	public boolean equals(SpecInterpretation t1, SpecInterpretation t2)
	{
		return (t1.getStart() == t2.getStart());
	}
}