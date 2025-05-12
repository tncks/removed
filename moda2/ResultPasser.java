package moda;

import modi.*;
import msutil.PGraph;
import processedDB.TagTrie;

import java.util.ArrayList;

public class ResultPasser {

    public ResultPasser() {
    }


    public ArrayList<AnsPeptide> dynamicMODeye(TagTrie dynamicDB, PGraph graph, TagPool tPool) {
        SpectrumAnalyzer szer = new SpectrumAnalyzer();
        MatchedTagPool matchedList = szer.extendedBuildMatchedTagPool(tPool, graph.getCorrectedMW(),
                dynamicDB, Mutables.protease, Constants.numberOfEnzymaticTermini);

        TagChainPool tcPool = new TagChainPool();
        tcPool.putAll(szer.buildTagChain(matchedList));
        tcPool.discardPoorTagChain();

        boolean specAnnotated = false;
        if (tcPool.size() != 0) {
            specAnnotated = szer.interpretTagChain(Mutables.variableModifications, tcPool, graph);
        }

        ArrayList<AnsPeptide> cands = new ArrayList<>();
        if (tcPool.size() != 0 && specAnnotated) {
            cands = tcPool.getAnswerPeptides(graph);
        }
        return cands;
    }


}
