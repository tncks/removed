package modi;

import java.util.LinkedList;

public class PeptideMatchInfo {
	int						startPos, endPos;
	Sequence				sequence;
	LinkedList<MatchedTag>	matchedTagList;
	
	public	PeptideMatchInfo(int startPos, int endPos, Sequence sequence, LinkedList<MatchedTag> tagList)
	{
		this.startPos = startPos;
		this.endPos = endPos;
		this.sequence = sequence;
		this.matchedTagList = tagList;
	}
	
	int			getStartPos()	{ return startPos; }
	int			getEndPos()		{ return endPos; }
	Sequence	getSequence()	{ return sequence; }
	LinkedList<MatchedTag>	getMatchedTagList()	{ return matchedTagList; }
	
	public String toString() 
	{
		StringBuffer output = new StringBuffer();
		output.append(startPos + "," + endPos + " " + sequence + " ");
		for(MatchedTag tag : matchedTagList)
			output.append(tag.getSequence() + " ");
		output.append("\n");
		
		return output.toString();
	}
}
