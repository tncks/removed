package modi;

import java.util.LinkedList;

import msutil.PGraph;
import msutil.ProtCutter;

import processedDB.RetrivedPeptideMap;
import processedDB.TagTrie;

public class SpectrumAnalyzer {
	
	private	SpectrumAnalyzer() {}
	public	static TagPool buildTagPool( Spectrum sourceSpec ) {
		if( sourceSpec == null ) return null;

		for(int i=0; i<sourceSpec.size(); i++){
			if( Math.abs(sourceSpec.get(i).getMass() - sourceSpec.getPrecursor() ) < 2 ) {
				sourceSpec.remove(i);
				i--;			
				continue;
			}
			sourceSpec.get(i).setIndex(i);
		}//temporal
	
		sourceSpec.normalizeIntensityLocally();
		
		int extra = ( sourceSpec.getCharge() > 2 && Constants.INSTRUMENT_TYPE != Constants.msms_type.QTOF )? 2 : 0; 		
		sourceSpec.peakSelection(Constants.selectionWindowSize, Constants.minNumOfPeaksInWindow+extra );	
		TagPool primitiveTags = sourceSpec.generateTags(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain, Constants.massToleranceForDenovo);
	
		return primitiveTags;
	}

	public	static TagChainPool buildTagChain(MatchedTagPool matchedTags)
	{
		TagChainPool tagChainPool = new TagChainPool();
		tagChainPool.buildTagChainPool(matchedTags);
		return tagChainPool;
	}
	
	public 	static void discoveryNovelPTM( ScanCap scap, PTMDB ptmDB, TagChainPool tcPool, PtmMap mapTable, PGraph pg )
	{
		TagChain bestTC = tcPool.getBestTagChain();
		
		if( bestTC == null || bestTC.getScore() < bestTC.tagCoverage*0.7 || 
				bestTC.matchedPeptide.size() < 8 ||
				(double)bestTC.tagCoverage/bestTC.matchedPeptide.size() < 0.5 )
			return;

		int novel= 0;
		for( SpecInterpretation t : bestTC ){
			if( t instanceof Gap ){
				Gap gap= (Gap)t;
				if( !Constants.fEqual(gap.getOffset(), 0) )
					novel ++;
			}
		}
		
		if( novel > 2 ) // current limit : two PTMs per peptide
			return;

		PTM[] tempPtm = new PTM[2];
		int pIndex = 0, gappedmatch= 0;
		for( SpecInterpretation t : bestTC ){
			if( t instanceof Gap ){
				Gap gap= (Gap)t;
				if( !Constants.fEqual(gap.getOffset(), 0) ){
					int site= gap.predictModifiedSite();
					// 0:N-term, Pos:site, -1:N-termMat, -2:C-termMat, -10:false
					if( site == -10 )
						continue;
					gappedmatch++;						
					if( site > -1 ){
						if( ptmDB.isNovelPtm( site, bestTC.matchedPeptide.get(site).getIndex(), gap.getOffset()) ){
							if( site == 0 )
								tempPtm[pIndex++]= new PTM(site,"","",gap.getOffset(), 0, bestTC.matchedPeptide.get(0), PTMPosition.ANY_N_TERM);
							else
								tempPtm[pIndex++]= new PTM(site,"","",gap.getOffset(), 0, bestTC.matchedPeptide.get(site), PTMPosition.ANYWHERE);
						}
					}
				}
			}
		}
		if ( novel == gappedmatch && pIndex == 1 ){
				mapTable.hitPtmFreqTable( scap, bestTC.matchedPeptide, tempPtm[0], pg);
		}
	}
	
	public 	static boolean interpretTagChain(PTMDB ptmDB, TagChainPool tcPool, PGraph graph )
	{				
		Spectrum sourceSpectrum = null;
		boolean specAnnotated = false;
		
		for ( LinkedList<TagChain> tagChainList : tcPool.values() )
		{
			for (int k=0;k<tagChainList.size(); k++)
			{
				TagChain tc = tagChainList.get(k);
				
				boolean allGapAnnotated = true;
				if(sourceSpectrum == null) {
					sourceSpectrum = tc.sourceSpectrum;
				}
				Peptide pep = tc.getMatchedPeptide();
				for ( SpecInterpretation si : tc ) {
					if (!(si instanceof Gap)) continue;
					Gap gap = (Gap)si;
					PTMSearchResult interpretation = ptmDB.searchPTM( pep.subSequence(gap.getStart(), gap.getEnd()+1 ), 
							gap.getOffset(), gap.getPosition() );

					if( !interpretation.isInterpreted() ) {
						gap.setInterpreted(false);
						allGapAnnotated = false;
						tc.setAllGapAnnotated(false);
						break;
					}
					else gap.setInterpreted(true);
					
					gap.setInterpretation(interpretation, graph);
				}
				
				if( allGapAnnotated ){
					tc.setAllGapAnnotated(true);
					specAnnotated = true;		
				}
				else{
					tagChainList.remove(k);
					k--;
				}
			}
		}
		return specAnnotated;
	}
	
	public static MatchedTagPool extendedBuildMatchedTagPool(TagPool primitiveTags, double motherMass, 
			TagTrie ixPDB, ProtCutter enzyme, int NTT)
	{
		if( primitiveTags == null || ixPDB == null )
			return null;
		
		double minDelta = (Constants.minModifiedMass < 0)? Constants.minModifiedMass - Constants.gapTolerance : - Constants.gapTolerance;
		double maxDelta = (Constants.maxModifiedMass > 0)? Constants.maxModifiedMass + Constants.gapTolerance : + Constants.gapTolerance;
		TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

		int realTag = 0;
		double orbMass= motherMass - Constants.H2O;		
		RetrivedPeptideMap searchResults= new RetrivedPeptideMap();
		for(Tag tag : longTags){
			
			RetrivedPeptideMap bRes= ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, tag.getBIonNtermOffset()-Constants.NTERM_FIX_MOD, tag, 
					tag.getBIonCtermOffset()-Constants.CTERM_FIX_MOD, IonDirection.B_DIRECTION, minDelta, maxDelta, Constants.gapTolerance);
			searchResults.combine(bRes);
			
			Tag reverseTag= tag.reverseTag();
			RetrivedPeptideMap yRes= ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, reverseTag.getYIonNtermOffset()-Constants.NTERM_FIX_MOD, reverseTag, 
					reverseTag.getYIonCtermOffset()-Constants.CTERM_FIX_MOD, IonDirection.Y_DIRECTION, minDelta, maxDelta, Constants.gapTolerance);
			searchResults.combine(yRes);
			realTag++;
			
			if( realTag > Constants.MAX_TAG_SIZE*2 ) break;
		}
		return searchResults.convertToMatchedTagPool( primitiveTags.extract(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain) );
	}

}





