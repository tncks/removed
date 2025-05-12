package modi;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;

public class TagPool extends ArrayList<Tag> {
	
	public void setTagScores(){
		for(Tag t : this){
			t.setTagScore();
		}
	}
	
	public TagPool extractQualifiedTag(int length, int size) {

		class XXXComparator implements Comparator<Tag> {
			public int compare(Tag t1, Tag t2) {
				
				String s1 = t1.sequence().toString();
				String s2 = t2.sequence().toString();
				
				if( s1.compareTo(s2) > 0 ) return 1;
				else if( s1.compareTo(s2) < 0 ) return -1;

				if( t1.getFirst().getMass() > t2.getFirst().getMass() ) return 1;
				else if( t1.getFirst().getMass() < t2.getFirst().getMass() ) return -1;

				return 0;
			}
		}

		this.sort(new XXXComparator());
		String prevSeq = "XXX";
		double prevPek = 0;
		TagPool refinedList = new TagPool();
		TagPool tpList = new TagPool();
		for(int i=0; i<this.size(); i++) {
			if( this.get(i).size() < length+1 ) continue;
			double curPek = this.get(i).getFirst().getMass();
			String curSeq = this.get(i).sequence().toString();
			
			if( Math.abs( prevPek-curPek ) < 4 && prevSeq.compareTo(curSeq) == 0 ){
				tpList.add( this.get(i) );
			}
			else {
				tpList.sort(Collections.reverseOrder(new TagComparator()));
				for(int k=0; k<tpList.size() && k<2; k++) refinedList.add(tpList.get(k));
				tpList = new TagPool();
				tpList.add( this.get(i) );
				prevPek = curPek;
				prevSeq = curSeq;
			}
		}

		tpList.sort(Collections.reverseOrder(new TagComparator()));
		for(int k=0; k<tpList.size()&&k<2; k++) refinedList.add(tpList.get(k));

		refinedList.sort(new TagComparator());

		TagPool extracted = new TagPool ();
		int selectedTagSize=0;	
		for(int i=refinedList.size()-1; i>=0 && selectedTagSize<=size; i--) {
			Tag curTag = refinedList.get(i);
			selectedTagSize++;
			if( (Constants.Leu_indistinguishable_Ile && curTag.sequence().contains(AminoAcid.getAminoAcid('I'))) ||
					(Constants.Lys_indistinguishable_Qln && curTag.sequence().contains(AminoAcid.getAminoAcid('K'))) )
				selectedTagSize--;			
			extracted.add( curTag );
		}	
		return extracted;
	}
	
	public TagPool extractAbove(int length)
	{
		TagPool extracted = new TagPool ();
		Iterator<Tag> it = this.iterator();
		Tag cur;
		while(it.hasNext())
			if(((cur = it.next()).sequence().size() >= length))
				extracted.add(cur);
		return extracted;
	}

	public TagPool extract(int start, int end)	// start : inclusive, end : exclusive
	{
		TagPool extracted = new TagPool ();
		Iterator<Tag> it = this.iterator();
		Tag cur;
		while(it.hasNext()) 
		{
			cur = it.next();
			if(cur.sequence().size() >= start && cur.sequence().size() < end)
				extracted.add(cur);
		}
		return extracted;
	}
	
	
//	public static TagPool extendTags(TagPool target, TagPool source)
//	{
//		TagPool result = new TagPool();
//		for(int i=0; i<target.size(); i++)
//			for(int j=0; j<source.size(); j++)
//				if(Tag.extendable(target.get(i), source.get(j)))
//					result.add(Tag.merge(target.get(i), source.get(j)));
//
//		return result;
//	}
	
	public String toString()
	{		
		StringBuilder ret = new StringBuilder();
		Iterator<Tag> it = this.iterator();
		while(it.hasNext())
			ret.append(it.next().toString()).append("\n");
			
		return ret.toString();
	}
	
}

