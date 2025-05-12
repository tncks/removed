package processedDB;

import java.util.ArrayList;
import java.util.Collections;

import moda.ThreadLocalMutables;
import modi.Constants;
import modi.Mutables;
import msutil.*;

public class ChainTagPeptide extends MODPeptide {
	

	final ArrayList<SequenceTag> mTags;

	public ChainTagPeptide(int s, int e, SequenceTag tag){
		super(s, e);
		mTags = new ArrayList<>();
		mTags.add(tag);
	}
	
	public ArrayList<SequenceTag> getMatchedTags( ){ return mTags; }

	public boolean extend( TagPeptide merged, TagPeptide toMerge, TagTrie trie ){ //multi-mod

		double	a = ( Constants.maxModifiedMass < 0 )?  (ThreadLocalMutables.get().gapTolerance) : Constants.maxModifiedMass+(ThreadLocalMutables.get().gapTolerance);
		double	b = ( Constants.minModifiedMass > 0 )? -(ThreadLocalMutables.get().gapTolerance) : Constants.minModifiedMass-(ThreadLocalMutables.get().gapTolerance);
		double	shiftWindow = (a - b);
		
		if( merged.pStart <= toMerge.pStart && toMerge.pStart < merged.pEnd )
		{
			double offset = MSMass.getPepMass( trie.getPeptide(merged.pStart, toMerge.pStart) );
			if( offset <= shiftWindow && mTags.size() < 100 ) {
				this.pLeft = Math.max( this.pLeft,  toMerge.pLeft );
				this.pRight= Math.min( this.pRight, toMerge.pRight );
				this.pEnd  = toMerge.pEnd;
				
				SequenceTag tag = toMerge.mTag;
				tag.shiftTag( getNOffset(merged, trie)+offset, (toMerge.pStart - this.pStart) );
				this.mTags.add(tag);	
				
				return true;
			}
		}
		return false;
	}
	
	public double getNOffset( TagPeptide xp, TagTrie trie ){
		return msutil.MSMass.getPepMass( trie.getPeptide(this.pStart, xp.pStart) );		
	}
	
	public void arrangeTags(){
		Collections.sort( mTags );
		int initSize = mTags.size();
		for(int i=0; i<initSize-1; i++){					
			for(int j=i+1; j<initSize; j++){
				if( mTags.get(i).getEnd()+1 < mTags.get(j).getStart() ) break;
				if( mTags.get(i).extend( mTags.get(j) ) ){
					mTags.remove(j);		
					j--;
					initSize--;
                }
			}
		}			
	}
}
