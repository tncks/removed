package modi;

import java.text.DecimalFormat;

import msutil.MSMass;


public class Constants {


    public static final String DECOY_LABEL = "dec_";

    public static final double UNIT_MASS = 1.;

    public static final double Electron = 0.000549;
    public static final double Hydrogen = 1.007825035;
    public static final double Oxygen = 15.99491463;
    public static final double Nitrogen = 14.003074;
    public static final double Proton = Hydrogen - Electron;
    public static final double H2O = Hydrogen * 2 + Oxygen;
    public static final double NH3 = Hydrogen * 3 + Nitrogen;
    public static final double IsotopeSpace = 1.00235;


    public static final double B_ION_OFFSET = Proton;
    public static final double Y_ION_OFFSET = H2O + Proton;
    public static final double A_ION_OFFSET = Oxygen + 12.;
    public static final double IMM_OFFSET = -A_ION_OFFSET + Proton;

    public static final double minPeptideMass = 300.;
    public static final double maxPeptideMass = 5000.;

    public static final double minNormIntensity = 0.00;

    public static final int maxTagPerPept = 12;
    public static final int maxTagChainPerPept = 30;
    public static final int maxInterpretationPerGap = 10;
    public static final double selectionWindowSize = 70;
    public static final String UNIMOD_FILE_NAME = "unimod.xml";
    public static final double[] rNorm = {6,
            2.928968, 1.928968, 1.428968, 1.095635, 0.845635,
            0.645635, 0.478968, 0.336111, 0.211111, 0.100000};

    public static int targetDecoy = 0;

    public static String firstSearchProgram = "";

    public static double NTERM_FIX_MOD = 0;
    public static double CTERM_FIX_MOD = 0;


    public static int numberOfEnzymaticTermini = 2;
    public static int missCleavages = 2;

    public static int minNoOfC13 = 0;

    public static int rangeForIsotopeIncrement = 0;

    public static double alkylatedToCys = 0;
    public static String alkylationMethod;


    public static double PPMTolerance = 0;


    public static double minModifiedMass = -0.5;
    public static double maxModifiedMass = 0.5;
    public static boolean canBeModifiedOnFixedAA = false;


    public static int minNumOfPeaksInWindow = 4;


    public static boolean Leu_indistinguishable_Ile = true;
    public static boolean Lys_indistinguishable_Qln = true;


    public static int maxPTMPerGap = 2;
    public static int maxPTMPerPeptide = 4;

    public static String PTM_FILE_NAME = "PTMDB.xml";


    public static String isobaricTag = "";
    public static String enrichedModification = "";


    public enum spectra_format {
        PKL,
        DTA,
        MGF,
        MS2,
        MZXML,
        ZIPDTA,
    }

    public enum msms_type {
        QTOF,
        TRAP,
    }


    public static String engine;
    public static String engineVersion;
    public static String runDate;
    public static String runUser = "anonymous";
    public static String runTitle;

    public static String SPECTRUM_LOCAL_PATH;
    public static String PROTEIN_DB_LOCAL_PATH;
    public static String INSTRUMENT_NAME = "TRAP";
    public static msms_type INSTRUMENT_TYPE = msms_type.TRAP;
    public static spectra_format SPECTRA_FILE_TYPE = spectra_format.MGF;

    public static int MSResolution = 0;
    public static int MSMSResolution = 0;

    public static double tagChainPruningRate = 0.5;
    public static int minTagLength = 3;
    public static int MAX_TAG_SIZE = 50;
    public static int minTagLengthPeptideShouldContain = 3;


    public static void adjustParameters() {
        if (INSTRUMENT_TYPE == msms_type.QTOF) { // TOF
            Mutables.massToleranceForDenovo = 0.2;
            Constants.minNumOfPeaksInWindow = 4;
            rNorm[0] = 6;
        } else {
            Mutables.massToleranceForDenovo = (MSMSResolution == 0) ? 0.3 : 0.03;
            Constants.minNumOfPeaksInWindow = 4;
            rNorm[0] = 6;
        }
        if (Mutables.massToleranceForDenovo > Mutables.fragmentTolerance / 2)
            Mutables.massToleranceForDenovo = Mutables.fragmentTolerance / 2;
        if (Mutables.fragmentTolerance < 0.1) MSMSResolution = 1;
        Constants.Lys_indistinguishable_Qln = MSMass.isIndistinguishableAA('K', 'Q');
        Constants.Leu_indistinguishable_Ile = MSMass.isIndistinguishableAA('L', 'I');

        if (Constants.canBeModifiedOnFixedAA) {
            double fixedOff = -20;
            if (!Mutables.fixedModifications.isEmpty()) {
                for (PTM p : Mutables.fixedModifications) {
                    fixedOff -= p.getMassDifference();
                }
                if (fixedOff < Constants.minModifiedMass) Constants.minModifiedMass = fixedOff;
            }
        }
        Constants.MAX_TAG_SIZE = 100;
        Constants.minTagLength = 2;
        Constants.minTagLengthPeptideShouldContain = 3;
        Constants.tagChainPruningRate = 0.4;
    }


    public static String getString(double value) {
        return new DecimalFormat("#.###").format(value);
    }

    public static int getMaxPTMOccurrence(int seqLength) {
        if (seqLength > 10) return 1;
        return Constants.maxPTMPerGap;
    }


    public static double PPMtoDalton(double mass, double ppm) {
        return mass / 1000000 * ppm;
    }

    public static int round(double a) {
        if (a > 0) return (int) (a + 0.5);
        else return (int) (a - 0.5);
    }
}
