package processedDB;

import modi.IonDirection;
import modi.Mutables;
import modi.Tag;
import modi.Constants;

public class TagMatchToPept implements Comparable<TagMatchToPept> {
	
	final Tag 			matchedTag;
	final double 			nGap;
    final double cGap;
	final int 			staSite;
    int endSite;
	final IonDirection 	ir;
	
	public TagMatchToPept(Tag tag, double n, double c, int start, int end, IonDirection x){
		matchedTag= tag;
		nGap= n;
		cGap= c;
		staSite= start;
		endSite= end;
		ir= x;
	}

	public boolean extendable(TagMatchToPept x){
		
		if( Mutables.fEqual(this.nGap, x.nGap) ){
            return x.staSite >= this.staSite && x.staSite <= this.endSite + 1;
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


