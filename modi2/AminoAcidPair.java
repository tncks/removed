package modi;

import java.util.ArrayList;
import java.util.Collections;

public class AminoAcidPair implements Comparable<AminoAcidPair>
{
	static 
	{
		AminoAcidPairList = new ArrayList<>();
		makeAAPairTable();
	}
	
	private final AminoAcid	firstAA;
	private final AminoAcid	secondAA;
	private final double		mass;
	
	public	AminoAcidPair(AminoAcid firstAA, AminoAcid secondAA, double mass)
	{
		this.firstAA = firstAA;
		this.secondAA = secondAA;
		this.mass = mass;
	}

	public	int	compareTo(AminoAcidPair tag)
	{
		if(this.mass > tag.mass)
			return 1;
		else if(this.mass == tag.mass)
			return 0;
		else
			return -1;
	}
	
	public	String toString()
	{
		return "" + firstAA + secondAA + " : " + mass;
	}
	
	private static final ArrayList<AminoAcidPair> AminoAcidPairList;

	private	static void	makeAAPairTable()
	{
		AminoAcid firstAA, secondAA;
		
		firstAA = AminoAcid.getAminoAcid('P');
		for(int j=0; j<'Z'-'A'; j++)
		{
			secondAA = AminoAcid.getAminoAcid((char)('A'+j));
			if(secondAA == null)
				continue;
			AminoAcidPairList.add(new AminoAcidPair(firstAA, secondAA, firstAA.getMonoMass()+secondAA.getMonoMass()));
		}

		Collections.sort(AminoAcidPairList);
	}
}
