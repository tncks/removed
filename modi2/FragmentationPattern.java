package modi;

public class FragmentationPattern {

	public static int		getClassification(Spectrum sourceSpectrum, Peptide matchedPeptide)
	{
		int classification = 0;

		// charge
		int charge = sourceSpectrum.getCharge();
		
		// terminal
		PeptTerminal terminal = PeptTerminal.NON_TRYPTIC;	
		AminoAcid peptTerm;
		if((peptTerm = matchedPeptide.getLast()) == AminoAcid.getAminoAcid('R'))
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
