package modi;

public class NeutralLoss {
	public static final double AMMONIA = Constants.UNIT_MASS*17;
	public static final double WATER = Constants.UNIT_MASS*18;
	public static final double CARBOXYL_OXIDE = Constants.UNIT_MASS*28;
	public static final double DIOXIDE_SULFIDE = Constants.UNIT_MASS*64;
	
	public static final double[] arrNeutralLoss =
	{
		AMMONIA,
		WATER,
		CARBOXYL_OXIDE,
		DIOXIDE_SULFIDE
	};
	
	public static boolean isNeutralLoss(double offset1, double offset2)
	{
		for(int i=0; i<arrNeutralLoss.length; i++)
		{
			if(Constants.fEqual(Math.abs(offset1-offset2), arrNeutralLoss[i]))
				return true;
		}
		
		return false;
	}
}
