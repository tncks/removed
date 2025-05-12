package msutil;

import modi.Constants;
import modi.Mutables;

import java.util.Arrays;

@SuppressWarnings("unused")
public class ProtCutter {

    private static final int N_TERM = 1;
    private static final int C_TERM = 2;
    private static final int BOTH_TERM = 3;

    private final String name;
    private String cleave = "";
    private String term = "";
    private final boolean specified;

    private String nTermCleave = "", cTermCleave = "";
    private final int[] cleaveMap = new int[26];

    public ProtCutter(String n, String c, String t) {
        name = n;
        cleave = c;
        term = t;
        specified = true;

        for (int i = 0; i < term.length(); i++) {
            if (term.charAt(i) == 'C') cTermCleave = c;
            if (term.charAt(i) == 'N') nTermCleave = c;
        }
        setCleaveMap(nTermCleave, cTermCleave, true);
    }

    public ProtCutter(String na, String nc, String cc, boolean isEnzyme) {
        name = na;
        nTermCleave = nc;
        cTermCleave = cc;
        specified = isEnzyme;
        setCleaveMap(nTermCleave, cTermCleave, isEnzyme);

    }

    public void setCleaveMap(String nc, String cc, boolean isEnzyme) {

        if (!isEnzyme) {
            Arrays.fill(cleaveMap, C_TERM);
            Constants.numberOfEnzymaticTermini = 0;
        } else {
            for (int i = 0; i < nTermCleave.length(); i++) {
                cleaveMap[nTermCleave.charAt(i) - 'A'] = N_TERM;
            }
            for (int i = 0; i < cTermCleave.length(); i++) {
                {
                    if (cleaveMap[cTermCleave.charAt(i) - 'A'] == N_TERM)
                        cleaveMap[cTermCleave.charAt(i) - 'A'] = BOTH_TERM;
                    else cleaveMap[cTermCleave.charAt(i) - 'A'] = C_TERM;
                }
            }
        }
    }

    public String getName() {
        return name;
    }

    public String getCleave() {
        return cleave;
    }

    public String getNtermCleave() {
        return nTermCleave;
    }

    public String getCtermCleave() {
        return cTermCleave;
    }

    public String getTerm() {
        return term;
    }

    public boolean isSpecified() {
        return specified;
    }

    public boolean isCleavage(int prev_aa, int next_aa) {
        char TERMINAL_RES = '-';
        if (prev_aa == TERMINAL_RES || next_aa == TERMINAL_RES) return true;
        if (cleaveMap[prev_aa - 'A'] == C_TERM || cleaveMap[prev_aa - 'A'] == BOTH_TERM) return true;
        return cleaveMap[next_aa - 'A'] == N_TERM || cleaveMap[next_aa - 'A'] == BOTH_TERM;
    }


    public static ProtCutter getCutter(String name) {
        for (ProtCutter e : enzymeTable)
            if (e.getName().equalsIgnoreCase(name))
                return e;
        return null;
    }

    private static final ProtCutter[] enzymeTable =
            {
                    new ProtCutter("Trypsin", "KR", "C"),
                    new ProtCutter("Arg-C", "R", "C"),
                    new ProtCutter("Glu-C", "E", "C"),
                    new ProtCutter("Asp-N", "BD", "N"),
                    new ProtCutter("Asp-N_ambic", "DE", "N"),
                    new ProtCutter("Chymotrypsin", "FYWL", "C"),
                    new ProtCutter("CNBr", "M", "C"),
                    new ProtCutter("Formic_acid", "D", "C"),
                    new ProtCutter("Lys-C", "K", "C"),
                    new ProtCutter("Lys-C/P", "K", "C"),
                    new ProtCutter("PepsinA", "FL", "C"),
                    new ProtCutter("Tryp-CNBr", "KRM", "C"),
                    new ProtCutter("TrypChymo", "FYWLKR", "C"),
                    new ProtCutter("Trypsin/P", "KR", "C"),
                    new ProtCutter("V8-DE", "BDEZ", "C"),
                    new ProtCutter("V8-E", "EZ", "C"),
                    new ProtCutter("CNBr+Trypsin", "KRM", "C"),
            };//*/

}
