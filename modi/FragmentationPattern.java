package modi;

public class FragmentationPattern {
	public static double	getMatchedPeakWeight(double relativeIntensity, Peak theoPeak, int classification)
	{
		double matchCredit = 0;
		PeakProperty property = theoPeak.getPeakProperty();
		if(property == PeakProperty.B_ION)
		{
			assert(theoPeak instanceof MatchedPeak);
			MatchedPeak mPeak = (MatchedPeak)theoPeak;
			AminoAcid NTermResidue = mPeak.getNTermResidue();
			AminoAcid CTermResidue = mPeak.getCTermResidue();
			switch(classification) 
			{
			case 1:
				return 0.2;
			case 2:
				return 0.25;
			case 3:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.1)
						return 0.5;
				}
				break;
			case 4:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.1)
						return 0.5;
				}
				break;
			case 5:
				if(NTermResidue == AminoAcid.getAminoAcid('D'))
				{
					if(relativeIntensity > 0.4)
						return 2;
				}
				else if(NTermResidue == AminoAcid.getAminoAcid('E'))
				{
					if(relativeIntensity > 0.2)
						return 1;
				}
				break;
			case 6:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.2)
						return 0.5;
				}
				break;
			case 7:
				if(CTermResidue == AminoAcid.getAminoAcid('I') ||
						CTermResidue == AminoAcid.getAminoAcid('K') ||
						CTermResidue == AminoAcid.getAminoAcid('L'))
				{
					if(relativeIntensity > 0.25)
						return 0.7;
				}
				break;
			case 8:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.2)
						return 0.5;
				}
				break;
			}
			return 0.25;
		}
		else if(property == PeakProperty.Y_MINUS_NH3_ION)
			return 0.2;
		else if(property == PeakProperty.Y_MINUS_H2O_ION)
			return 0.2;
		else if(property == PeakProperty.Y_ION)
		{
			assert(theoPeak instanceof MatchedPeak);
			MatchedPeak mPeak = (MatchedPeak)theoPeak;
			AminoAcid NTermResidue = mPeak.getNTermResidue();
			AminoAcid CTermResidue = mPeak.getCTermResidue();
			switch(classification) 
			{
			case 1:	
				if(NTermResidue == AminoAcid.getAminoAcid('D'))
				{
					if(relativeIntensity > 0.5)
						return 3;
				}
				else if(NTermResidue == AminoAcid.getAminoAcid('E'))
				{
					if(relativeIntensity > 0.3)
						return 2;
				}
				break;
			case 2:
				if(NTermResidue == AminoAcid.getAminoAcid('D'))
				{
					if(relativeIntensity > 0.3)
						return 2;
				}
				else if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.5)
						return 3;
				}
				break;
			case 3:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.4)
						return 2.5;
				}
				break;
			case 4:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.4)
						return 2;
				}
				break;
			case 5:
				if(NTermResidue == AminoAcid.getAminoAcid('D'))
				{
					if(relativeIntensity > 0.4)
						return 2;
				}
				else if(NTermResidue == AminoAcid.getAminoAcid('E'))
				{
					if(relativeIntensity > 0.2)
						return 1.5;
				}
				break;
			case 6:
				if(CTermResidue == AminoAcid.getAminoAcid('P'))
				{
					if(relativeIntensity > 0.4)
						return 2.5;
				}
				break;
			case 7:
				if(NTermResidue == AminoAcid.getAminoAcid('D'))
				{
					if(relativeIntensity > 0.4)
						return 2.5;
				}
				if(NTermResidue == AminoAcid.getAminoAcid('E'))
				{
					if(relativeIntensity > 0.25)
						return 2;
				}
				break;
			case 8:
				break;
			}
			return 1;
		}
		else if(property == PeakProperty.A_ION)
			return 0.1;
		else
			return 0;
	}
	
	private static double	getYIonIntensityLevel(double relativeIntensity, int classification, 
			AminoAcid NTermResidue, AminoAcid CTermResidue)
	{
		double matchWeight = 0;
		switch(classification)
		{
		case 1:
		}
		return matchWeight;
	}
	public static int		getClassification(Spectrum sourceSpectrum, Peptide matchedPeptide)
	{
		int classification = 0;

		// charge
		int charge = sourceSpectrum.getCharge();
		
		// terminal
		PeptTerminal terminal = PeptTerminal.NON_TRYPTIC;	
		AminoAcid peptTerm;
		if((peptTerm = matchedPeptide.get(matchedPeptide.size()-1)) == AminoAcid.getAminoAcid('R'))
			terminal = PeptTerminal.ARG;
		else if(peptTerm == AminoAcid.getAminoAcid('K'))
			terminal = PeptTerminal.ARG;
		

		// basic residues
		int numOfH = 0, numOfK = 0, numOfR = 0, numOfP = 0, numOfDE = 0;
		boolean KRprecedesP = false, KRprecedesPprecedesDE = false;
		for(AminoAcid aa : matchedPeptide)
		{
			if(aa == AminoAcid.getAminoAcid('H'))
				numOfH++;
			else if(aa == AminoAcid.getAminoAcid('K'))
				numOfK++;
			else if(aa == AminoAcid.getAminoAcid('R'))
				numOfR++;
			else if(aa == AminoAcid.getAminoAcid('P'))
			{
				if(numOfK > 0 || numOfR > 0)
					KRprecedesP = true;
				numOfP++;
			}
			else if(aa == AminoAcid.getAminoAcid('D') || aa == AminoAcid.getAminoAcid('E'))
			{
				if(KRprecedesP == true)
					KRprecedesPprecedesDE = true;
				numOfDE++;
			}
		}
		
		if(numOfP == 1 && numOfH == 0)
		{
			if(charge == 1)
			{
				if(terminal == PeptTerminal.ARG)
					classification = 1;
				else if(terminal == PeptTerminal.LYS)
					classification = 2;
			}
			else if(charge == 2)
			{
				if(terminal == PeptTerminal.ARG)
					classification = 3;
				else if(terminal == PeptTerminal.LYS)
					classification = 4;
			}
		}
		else if(KRprecedesPprecedesDE && numOfH == 0 && numOfDE == 1 && charge == 2)
		{
			if(terminal == PeptTerminal.ARG && numOfR == 2)
				classification = 5;
			else if(terminal == PeptTerminal.LYS && numOfK == 2)
				classification = 6;
		}
		else if(numOfP == 0 && numOfH == 0 && charge == 1 && terminal == PeptTerminal.LYS)
			classification = 7;
		else if(numOfP == 0 && numOfH == 0 && charge == 2 && (terminal == PeptTerminal.LYS || terminal == PeptTerminal.ARG))
			classification = 8;
		else
			classification = 9;
		
		return classification;
	}
}
