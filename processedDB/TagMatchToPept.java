package processedDB;

import modi.IonDirection;
import modi.Tag;
import modi.Constants;

public class TagMatchToPept implements Comparable<TagMatchToPept> {
	
	Tag 			matchedTag;
	double 			nGap, cGap;
	int 			staSite, endSite;	
	IonDirection 	ir;
	
	public TagMatchToPept(Tag tag, double n, double c, int start, int end, IonDirection x){
		matchedTag= tag;
		nGap= n;
		cGap= c;
		staSite= start;
		endSite= end;
		ir= x;
	}
	
	public Tag getMatchedTag(){ return matchedTag; }
	public double getDelta(){ return nGap+cGap; }
	
	public boolean extendable(TagMatchToPept x){
		
		if( Constants.fEqual(this.nGap, x.nGap) ){
			if( x.staSite >= this.staSite && x.staSite<= this.endSite+1 )
				return true;
		}		
		return false;
	}
	
	public boolean extend(TagMatchToPept x) {
		if( extendable(x) ){
			this.endSite= Math.max( this.endSite, x.endSite );
			return true;
		}	
		return false;
	}
	
	public int compareTo(TagMatchToPept x){	// default comparator : mass
		if( x.staSite < this.staSite ) return 1;
		else if( x.staSite > this.staSite ) return -1;
		else return 0;
	}	
	
}


