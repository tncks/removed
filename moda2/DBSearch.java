package moda;

import java.util.LinkedList;

import modi.AminoAcid;
import modi.Constants;
import modi.Mutables;
import modi.PeakProperty;
import modi.Tag;
import modi.TagPool;
import processedDB.CandidateContainer;
import processedDB.MODPeptide;
import processedDB.TagPeptide;
import processedDB.TagTrie;

public class DBSearch {
	public	static CandidateContainer construct_onemod_cpool(TagPool primitiveTags, double motherMass, TagTrie ixPDB) {

		if( primitiveTags == null || ixPDB == null )
			return null;

		Mutables mut = ThreadLocalMutables.get();
		double minDelta = (Constants.minModifiedMass < 0)? Constants.minModifiedMass - mut.gapTolerance : - mut.gapTolerance;
		double maxDelta = (Constants.maxModifiedMass > 0)? Constants.maxModifiedMass + mut.gapTolerance : + mut.gapTolerance;

		TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

		LinkedList<MODPeptide> cpool= new LinkedList<MODPeptide>();

		int realTag = 0, redunTag = 0;
		for(Tag tag : longTags){
			if( tag.get(0).getPeakProperty() != PeakProperty.C_TERM_Y_ION_ONLY &&
					tag.get(3).getPeakProperty() != PeakProperty.N_TERM_Y_ION_ONLY )
			{
				LinkedList<MODPeptide> bRes= ixPDB.getOneModPeptides(tag.getBIonNtermOffset()-Constants.NTERM_FIX_MOD, tag.sequence().toString(),
						tag.getBIonCtermOffset()-Constants.CTERM_FIX_MOD, 0, minDelta, maxDelta, mut.gapTolerance);
				cpool.addAll(bRes);
			}
			if( tag.get(0).getPeakProperty() != PeakProperty.N_TERM_B_ION_ONLY &&
					tag.get(3).getPeakProperty() != PeakProperty.C_TERM_B_ION_ONLY )
			{
				Tag reverseTag= tag.reverseTag();
				LinkedList<MODPeptide> yRes= ixPDB.getOneModPeptides(reverseTag.getYIonNtermOffset()-Constants.NTERM_FIX_MOD, reverseTag.sequence().toString(),
						reverseTag.getYIonCtermOffset()-Constants.CTERM_FIX_MOD, 1, minDelta, maxDelta, mut.gapTolerance);
				cpool.addAll(yRes);
			}
			realTag++;
			redunTag++;

			if( tag.sequence().contains(AminoAcid.getAminoAcid('I')) || tag.sequence().contains(AminoAcid.getAminoAcid('K')) )
				redunTag--;
			if( realTag > MODaConst.maxTagPoolSize*2 || redunTag > MODaConst.maxTagPoolSize ) break;//*/
		}

		return new CandidateContainer(cpool);
	}

	public	static CandidateContainer construct_multimod_cpool(TagPool primitiveTags, double motherMass, TagTrie ixPDB) {

		if( primitiveTags == null || ixPDB == null )
			return null;

		Mutables mut = ThreadLocalMutables.get();
		double minDelta = (Constants.minModifiedMass < 0)? Constants.minModifiedMass - mut.gapTolerance : - mut.gapTolerance;
		double maxDelta = (Constants.maxModifiedMass > 0)? Constants.maxModifiedMass + mut.gapTolerance : + mut.gapTolerance;

		TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

		LinkedList<TagPeptide> cpool= new LinkedList<TagPeptide>();
		int realTag = 0, redunTag = 0;
		for(Tag tag : longTags){
			if( tag.get(0).getPeakProperty() != PeakProperty.C_TERM_Y_ION_ONLY &&
					tag.get(3).getPeakProperty() != PeakProperty.N_TERM_Y_ION_ONLY )
			{
				LinkedList<TagPeptide> bRes= ixPDB.getMultiModPeptides(tag.getBIonNtermOffset()-Constants.NTERM_FIX_MOD, tag.sequence().toString(),
						tag.getBIonCtermOffset()-Constants.CTERM_FIX_MOD, 0, minDelta, maxDelta, mut.gapTolerance);
				cpool.addAll(bRes);
			}
			if( tag.get(0).getPeakProperty() != PeakProperty.N_TERM_B_ION_ONLY &&
					tag.get(3).getPeakProperty() != PeakProperty.C_TERM_B_ION_ONLY )
			{
				Tag reverseTag= tag.reverseTag();
				LinkedList<TagPeptide> yRes= ixPDB.getMultiModPeptides(reverseTag.getYIonNtermOffset()-Constants.NTERM_FIX_MOD, reverseTag.sequence().toString(),
						reverseTag.getYIonCtermOffset()-Constants.CTERM_FIX_MOD, 1, minDelta, maxDelta, mut.gapTolerance);
				cpool.addAll(yRes);
			}
			realTag++;
			redunTag++;

			if( tag.sequence().contains(AminoAcid.getAminoAcid('I')) || tag.sequence().contains(AminoAcid.getAminoAcid('K')) )
				redunTag--;
			if( realTag > MODaConst.maxTagPoolSize*2 || redunTag > MODaConst.maxTagPoolSize ) break;
		}

		return new CandidateContainer(cpool, ixPDB);
	}


//
//	public	static CandidateContainer construct_onemod_cpool(TagPool primitiveTags, TagTrie ixPDB) {
//
//		if( primitiveTags == null || ixPDB == null )
//			return null;
//
//		double minDelta = (Mutables.minModifiedMass < 0)? Mutables.minModifiedMass - Mutables.gapTolerance : - Mutables.gapTolerance;
//		double maxDelta = (Mutables.maxModifiedMass > 0)? Mutables.maxModifiedMass + Mutables.gapTolerance : + Mutables.gapTolerance;
//
//		TagPool longTags = primitiveTags.extractAbove(Mutables.minTagLengthPeptideShouldContain);
//
//		LinkedList<MODPeptide> cpool= new LinkedList<>();
//
//		int realTag = 0, redunTag = 0;
//		for(Tag tag : longTags){
//			if( tag.get(0).getPeakProperty() != PeakProperty.C_TERM_Y_ION_ONLY &&
//					tag.get(3).getPeakProperty() != PeakProperty.N_TERM_Y_ION_ONLY )
//			{
//				LinkedList<MODPeptide> bRes= ixPDB.getOneModPeptides(tag.getBIonNtermOffset()-Mutables.NTERM_FIX_MOD, tag.sequence().toString(),
//						tag.getBIonCtermOffset()-Mutables.CTERM_FIX_MOD, 0, minDelta, maxDelta, Mutables.gapTolerance);
//				cpool.addAll(bRes);
//			}
//			if( tag.get(0).getPeakProperty() != PeakProperty.N_TERM_B_ION_ONLY &&
//					tag.get(3).getPeakProperty() != PeakProperty.C_TERM_B_ION_ONLY )
//			{
//				Tag reverseTag= tag.reverseTag();
//				LinkedList<MODPeptide> yRes= ixPDB.getOneModPeptides(reverseTag.getYIonNtermOffset()-Mutables.NTERM_FIX_MOD, reverseTag.sequence().toString(),
//						reverseTag.getYIonCtermOffset()-Mutables.CTERM_FIX_MOD, 1, minDelta, maxDelta, Mutables.gapTolerance);
//				cpool.addAll(yRes);
//			}
//			realTag++;
//			redunTag++;
//
//			if( tag.sequence().contains(AminoAcid.getAminoAcid('I')) || tag.sequence().contains(AminoAcid.getAminoAcid('K')) )
//				redunTag--;
//			if( realTag > MODaConst.maxTagPoolSize*2 || redunTag > MODaConst.maxTagPoolSize ) break;//*/
//		}
//
//		return new CandidateContainer(cpool);
//	}
//
//	public	static CandidateContainer construct_multimod_cpool(TagPool primitiveTags, TagTrie ixPDB) {
//
//		if( primitiveTags == null || ixPDB == null )
//			return null;
//
//		double minDelta = (Mutables.minModifiedMass < 0)? Mutables.minModifiedMass - Mutables.gapTolerance : - Mutables.gapTolerance;
//		double maxDelta = (Mutables.maxModifiedMass > 0)? Mutables.maxModifiedMass + Mutables.gapTolerance : + Mutables.gapTolerance;
//
//		TagPool longTags = primitiveTags.extractAbove(Mutables.minTagLengthPeptideShouldContain);
//
//		LinkedList<TagPeptide> cpool= new LinkedList<>();
//		int realTag = 0, redunTag = 0;
//		for(Tag tag : longTags){
//			if( tag.get(0).getPeakProperty() != PeakProperty.C_TERM_Y_ION_ONLY &&
//					tag.get(3).getPeakProperty() != PeakProperty.N_TERM_Y_ION_ONLY )
//			{
//				LinkedList<TagPeptide> bRes= ixPDB.getMultiModPeptides(tag.getBIonNtermOffset()-Mutables.NTERM_FIX_MOD, tag.sequence().toString(),
//						tag.getBIonCtermOffset()-Mutables.CTERM_FIX_MOD, 0, minDelta, maxDelta, Mutables.gapTolerance);
//				cpool.addAll(bRes);
//			}
//			if( tag.get(0).getPeakProperty() != PeakProperty.N_TERM_B_ION_ONLY &&
//					tag.get(3).getPeakProperty() != PeakProperty.C_TERM_B_ION_ONLY )
//			{
//				Tag reverseTag= tag.reverseTag();
//				LinkedList<TagPeptide> yRes= ixPDB.getMultiModPeptides(reverseTag.getYIonNtermOffset()-Mutables.NTERM_FIX_MOD, reverseTag.sequence().toString(),
//						reverseTag.getYIonCtermOffset()-Mutables.CTERM_FIX_MOD, 1, minDelta, maxDelta, Mutables.gapTolerance);
//				cpool.addAll(yRes);
//			}
//			realTag++;
//			redunTag++;
//
//			if( tag.sequence().contains(AminoAcid.getAminoAcid('I')) || tag.sequence().contains(AminoAcid.getAminoAcid('K')) )
//				redunTag--;
//			if( realTag > MODaConst.maxTagPoolSize*2 || redunTag > MODaConst.maxTagPoolSize ) break;
//		}
//
//		return new CandidateContainer(cpool, ixPDB);
//	}
//
}


