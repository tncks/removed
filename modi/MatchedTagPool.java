package modi;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

public class MatchedTagPool extends TreeMap<Peptide, LinkedList<MatchedTag>>{
	public int add(MatchedTag tag)
	{
		if( !tag.checkCompatibility() )
			return 0;
		// New peptide tag

		Peptide matchedPeptide = tag.getMatchedPeptide();
		LinkedList<MatchedTag> tagList = this.get( matchedPeptide );
		if( tagList == null )
		{
			tagList = new LinkedList<MatchedTag>();
			tagList.add(tag);
			put( tag.getMatchedPeptide(), tagList );
			return 1;
		}

		tagList.add(tag);
		return tagList.size();
	}
	
	public	void addShortTags(TagPool primitiveTags, double motherMass, DoublePair motherMassTolerance, DoublePair offsetTolerance)
	{
		// match short tags to matched database
		PeptideDBForPTMSearchi matchedDB = new PeptideDBForPTMSearchi();
		try {
			matchedDB.construct(this.keySet(), Constants.minTagLength, Constants.minTagLengthPeptideShouldContain-1, 0, 0, 0);
		} catch(Exception e)
		{
			e.printStackTrace();
		}
		
		TagPool shortTags = primitiveTags.extract(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain);
		
		for(Tag tag : shortTags)
		{
			if(!tag.isYOnly())			// b direction
			{
				ArrayList<PeptideDBHit> bResults = matchedDB.search(
						motherMass, 
						tag.sequence(), 
						tag.getBIonNtermOffset(), 
						tag.getBIonCtermOffset(),
						motherMassTolerance,
						offsetTolerance
						);
				
				for(PeptideDBHit hit : bResults)
					if(!(tag.isCtermOnly() && !hit.isCtermIon() || tag.isNtermOnly() && !hit.isNtermIon()))
					{
						MatchedTag newTag = new MatchedTag(tag, hit, IonDirection.B_DIRECTION);
						this.add(newTag);
					}
			}
			if(!tag.isBOnly())			// y direction
			{
				Tag reverseTag = tag.reverseTag();
				ArrayList<PeptideDBHit> yResults = matchedDB.search(
						motherMass, 
						tag.reverseSequence(), 
						tag.getYIonNtermOffset(), 
						tag.getYIonCtermOffset(),
						motherMassTolerance,
						offsetTolerance
						);
				for(PeptideDBHit hit : yResults)
					if(!(tag.isCtermOnly() && !hit.isCtermIon() || tag.isNtermOnly() && !hit.isNtermIon()))
					{
						MatchedTag newTag = new MatchedTag(tag, hit, IonDirection.Y_DIRECTION);
						this.add(newTag);
					}
			}
		}
	}//*/
	
	public String toString() 
	{
		StringBuffer output = new StringBuffer();
		output.append("Matched Tag Pool");
		Iterator<Map.Entry<Peptide, LinkedList<MatchedTag>>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<MatchedTag>> entry = it.next();
			output.append(entry.getKey()).append("\n");
			output.append(entry.getValue()).append("\n");
		}
		return output.toString();
	}
	
}
