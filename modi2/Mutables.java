package modi;

import msutil.ProtCutter;



public class Mutables {
    public static ProtCutter protease = ProtCutter.getCutter("Trypsin");
    public static double fragmentTolerance = 0.6; // never change after set param invocation
    public static PTMDB variableModifications;
    public static PTMDB fixedModifications;
    public static double[] reporterMassOfIsobaricTag = null;
    public static double massToleranceForDenovo = 0.3;

    // Configuration values that should be copied to each thread
    public static int    DEFAULT_MAXNOOFC13 = 0;
    public static double DEFAULT_PRECURSORTOLERANCE = 0.5;
    public static double DEFAULT_PRECURSORACCURACY = 0.5;
    public static double DEFAULT_GAPTOLERANCE = 0.6;
    public static double DEFAULT_GAPACCURACY = 1.6;
    public static double DEFAULT_NONMODIFIEDDELTA = massToleranceForDenovo;

    // Instance fields for thread-local storage
    public int maxNoOfC13 = DEFAULT_MAXNOOFC13;
    public double precursorTolerance = DEFAULT_PRECURSORTOLERANCE;
    public double precursorAccuracy = DEFAULT_PRECURSORACCURACY;
    public double gapTolerance = DEFAULT_GAPTOLERANCE;
    public double gapAccuracy = DEFAULT_GAPACCURACY;
    public double nonModifiedDelta = DEFAULT_NONMODIFIEDDELTA;
    // TODO: Getters for the instance fields



    public Mutables() {
    }




    public static boolean fEqual(double v1, double v2) {
        return Math.abs(v1 - v2) <= Mutables.fragmentTolerance;
    }






}
