package moda;

import modi.Constants;
import modi.Mutables;
import msutil.MSMass;

public class MODaConst {

	public static final int minPeaksCount = 4;
	
	static final int maxTagPoolSize = 50;
	
	static final double 	baseScore = -10000.;
	static final double 	isotopeUnit = ( Constants.MSMSResolution==1 )? Constants.IsotopeSpace : 1.0;
	static final PtmMassUnit ptmUnit = ( Constants.MSMSResolution==1 )? new PtmMassUnitForHighRES() : new PtmMassUnit();
	
	static final double 	maxIsotopeError =  -1*isotopeUnit;
	static final int	 	maxIntIsotopeError =  -1;
	static final int 		isotopePointsToBeCorrected = ( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF )? 2 : 3;
	
	static final double 	TERMINALMOD = Constants.NTERM_FIX_MOD + Constants.CTERM_FIX_MOD;
	static final double 	minimumDistance = MSMass.getMinAAMass() - Mutables.fragmentTolerance;
	 
	static final double 	minCutOfProb = 0.9;
	static final int		minConsecutiveAALength = 5;
}
