package modi;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.TreeMap;

import moda.ThreadLocalMutables;
import msutil.*;

public class TagChainPool extends TreeMap<Peptide, LinkedList<TagChain>> {

	public void buildTagChainPool(MatchedTagPool matchedTags)
	{
		if( matchedTags == null || matchedTags.isEmpty()){
			return;
		}
		
		LinkedList<TagChain> tcList;
		
		Iterator<Map.Entry<Peptide, LinkedList<MatchedTag>>> it = matchedTags.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<MatchedTag>> entry = it.next();
			tcList = TagChain.buildTagChainList(entry);
			if(tcList.isEmpty()){ continue; }
			Peptide matchedPeptide = new Peptide(entry.getKey());
			put(matchedPeptide, tcList);
		}
	}

	public void discardPoorTagChain()	// should be modified
	{
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();

		// calculating best score
		double bestScore = 0;
		// TagChain bet = null;
		while( it.hasNext() ) {
			for(TagChain tc : it.next().getValue())
			{				
				double curScore = tc.getScore();
				if( bestScore < curScore ){
					bestScore = curScore;
					// bet = tc;
				}
			}
		}
	
		// discarding tag chain		
		it = this.entrySet().iterator();
		ArrayList<Peptide> deletePepList = new ArrayList<>();
		
		while(it.hasNext())
		{
			LinkedList<TagChain> tcList = it.next().getValue();
			Peptide pep = tcList.getFirst().getMatchedPeptide();
			ListIterator<TagChain> listIt = tcList.listIterator();
			while(listIt.hasNext()) {
				TagChain tc = listIt.next();
				if( tc.getScore() <= bestScore * Constants.tagChainPruningRate ){
					listIt.remove();
				}	
			}
			
			if(tcList.isEmpty())
				deletePepList.add(pep);
		}
		
		for( Peptide pep : deletePepList )
			this.remove(pep);	
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		output.append("Matched Tag Pool");
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<TagChain>> entry = it.next();
			output.append("\n").append(entry.getKey()).append("\n");
			for(TagChain tc : entry.getValue())
			{
				output.append(tc).append("\n");
			}
		}
		return output.toString();
	}

	boolean isWithinTolerance(double calc, double obsv, double tol) {

		if (Constants.minNoOfC13 == 0 && (ThreadLocalMutables.get().maxNoOfC13) == 0) {
			return !(Math.abs(calc - obsv) > tol);
		} else {
			double tempError = obsv - calc;
			int isoerr = Constants.round(tempError / Constants.IsotopeSpace);
			if (isoerr < Constants.minNoOfC13 || (ThreadLocalMutables.get().maxNoOfC13) < isoerr) return false;
			return !(Math.abs(tempError - isoerr * Constants.IsotopeSpace) > (ThreadLocalMutables.get().precursorAccuracy));
		}
	}

	int getModEyeRankScore( String peptide, double[] ptms, PGraph graph ){
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, graph);
		else iGraph= new TRAPGraph(peptide, ptms, graph);

		if( !isWithinTolerance(iGraph.getCalculatedMW(), graph.getObservedMW(), (ThreadLocalMutables.get().precursorTolerance)) )
			return -1;

		iGraph.setScore(graph);
		return iGraph.getRankScore();
	}
	
	public ArrayList<AnsPeptide> getAnswerPeptides( PGraph graph ){		
		AnsHeap answerPepts = new AnsHeap();
		
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> entries = this.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Peptide, LinkedList<TagChain>> entry = entries.next();
			String pept = entry.getKey().toString();
			LinkedList<TagChain> alignedTagChainList = entry.getValue();
			if(alignedTagChainList.isEmpty()) continue;
			
			HashSet<PTMCombination> ptmComb = new HashSet<>();
			for( TagChain tc : alignedTagChainList ){
				if( tc.allGapAnnotated ) ptmComb.addAll( tc.getPTMCombination() );
			}
			
			Iterator<PTMCombination> iter = ptmComb.iterator();
			while( iter.hasNext() ){
				PTMCombination p = iter.next();
	
				int s = getModEyeRankScore(pept, p.ptms, graph);
				if( s < 0 ) continue;
				AnsPeptide candidate = new AnsPeptide(entry.getKey(), p.ptmComb, p.ptms, p.ptmList, s);
				answerPepts.add( candidate ); 
			}
		}
		
		return answerPepts.getFinalList(graph);
	}
	
}
	
class TagChainListComparator implements Comparator<Map.Entry<Peptide, LinkedList<TagChain>>>
{
	public int compare(Map.Entry<Peptide, LinkedList<TagChain>> tc1, Map.Entry<Peptide, LinkedList<TagChain>> tc2)
	{
		LinkedList<TagChain> tcList1 = tc1.getValue();
		LinkedList<TagChain> tcList2 = tc2.getValue();
		assert(tcList1 != null && tcList2 != null); 
		if(tcList1.isEmpty() && tcList2.isEmpty())
			return 0;
		else if(tcList1.isEmpty())
			return 1;
		else if(tcList2.isEmpty())
			return -1;
		
		if(tcList1 == null || tcList2 == null || tcList1.isEmpty() || tcList2.isEmpty())
			return 0;
		double score1 = tcList1.getFirst().getScore();
		double score2 = tcList2.getFirst().getScore();
		if(score1 > score2)
			return -1;
		else if(score1 == score2)
			return 0;
		else
			return 1;
	}
	
	public boolean equal(LinkedList<TagChain> tcList1, LinkedList<TagChain> tcList2)
	{
		assert(tcList1 != null && tcList2 != null && !tcList1.isEmpty() && !tcList2.isEmpty());
		return tcList1.getFirst().getScore() == tcList2.getFirst().getScore();
	}
	
}
