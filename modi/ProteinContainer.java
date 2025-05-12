package modi;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;


public class ProteinContainer extends ArrayList<Protein> {
	int id = 0; // Temporary name
	int sizeofResidues=0;

	public int getSizeOfResidues(){ return sizeofResidues; }
	public void parseFromFile( String name ) throws IOException
	{
		File file = new File(name);
		if( file.isDirectory() )
			parseFromFastaDir(name);
		else parseFromFasta(name);
	}
	public void parseFromFastaDir( String dirName ) throws IOException
	{

		File protDir = new File(dirName);
		
		if(!protDir.isDirectory())
			return;

		for(File file : protDir.listFiles())
		{
			if(file.getPath().contains(".fasta"))
				parseFromFasta(file.getPath());
		}
	}
	
	public int parseFromFasta( String fileName ) throws IOException
	{

		BufferedReader in;
		try {
			in = new BufferedReader( new FileReader(fileName) );
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return -1;
		}
		
		String s = in.readLine();
		while (true)
		{
			while (s!=null && s.length()==0)
			{
				s=in.readLine();
			}
			
			if (s==null) break;
			
			if (s.charAt(0)!='>')
				return -1;				

			Protein protein = new Protein(s, id++);
			
			while (true)
			{
				s = in.readLine();
				if (s==null || s.length()==0 || s.charAt(0)=='>') break;
				
				for (int i=0; i<s.length(); i++)
				{
					AminoAcid aa = AminoAcid.getAminoAcid(s.charAt(i));					
					if (aa==null){
						// System.out.printf("(Warning)Invalid AminoAcid residue found : %c\n", s.charAt(i));
						// return -1;
					}
					else{
						protein.add(aa);
					}
				}
			}
			this.add(protein);
			sizeofResidues += protein.size();
		}
		
		System.out.println( "[Protein DB] Entries: " + this.size() +" Residues: " + sizeofResidues);
		
		in.close();
		return 0;
	}
	public void writeToFastaFile(File fileName)
	{
		try {
			PrintWriter writer = new PrintWriter(new FileOutputStream(fileName));
			for(Protein prot : this)
			{
				writer.println(prot.getAnnotaion());
				writer.println(prot.getSequence());
			}
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public char[] getTeminalResiduesOfPeptide(int id, int s, int t){
		char[] nc = new char[2];
		if( s == 0 ) nc[0] = '-';
		else  nc[0]= this.get(id).get(s-1).getResidue();
		
		if( t == this.get(id).size()-1 ) nc[1] = '-';
		else nc[1]= this.get(id).get(t+1).getResidue();
		
		return nc;
	}
	
	public char[] getTeminalResiduesOfPeptide(Peptide pept){
		int id = pept.getSrcProteinInfo().get(0).getSrcProteinID();
		int s = pept.getSrcProteinInfo().get(0).getStartPos();
		int t = pept.getSrcProteinInfo().get(0).getEndPos();
		char[] nc = new char[2];
		if( s == 0 ) nc[0] = '-';
		else  nc[0]= this.get(id).get(s-1).getResidue();
		
		if( t == this.get(id).size()-1 ) nc[1] = '-';
		else nc[1]= this.get(id).get(t+1).getResidue();
		
		return nc;
	}

	public String getProteinIdentity(Peptide pept){
		return getProteinIdentity( pept.getSrcProteinInfo().get(0) );	
	}
	public String getProteinIdentity(SrcProteinInfo pept){
		return this.get(pept.getSrcProteinID()).annotation;	
	}
	public String getCuratedPeptide(Peptide pept){
		char[] nc = getTeminalResiduesOfPeptide(pept);	
		return String.format("%c.%s.%c", nc[0], pept.toString(), nc[1]);
	}
}