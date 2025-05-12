package modi;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

@SuppressWarnings("unused")
public class PeptideDBForPTMSearchi extends HashMap<Sequence, Peptide> {
	private final TreeNode		root = new TreeNode();
	protected double		minMass;
	protected double		maxMass;
	protected int			miscleavage;
	protected int			maxQuerySeqLength;
	protected int			minQuerySeqLength;
	
	// Constant values(might be changed to read from setting file)
	private static final int		minPeptideSeqLength		= 4;
	private static final int		maxPeptideSeqLength		= 50;
	
	private static class TreeNode implements Serializable
	{
		final TreeNode[]				childs = new TreeNode[AminoAcid.getIndexSize()];
		final ArrayList<PeptideDBHit>	bucket = new ArrayList<>();
		
		// Get hit results recursively(including childs)
		public void xxx_getHitResult(double motherMass, double nTermOffset, double cTermOffset, DoublePair motherMassTolerance, DoublePair offsetTolerance, ArrayList<PeptideDBHit> results)
		{
			// check if offset fits the max offset difference, and add them to result
			
			for (PeptideDBHit dbHit : bucket)
			{
				double pepMass = dbHit.getHitPeptide().getMonoMass()+Constants.H2O;
				if ( motherMassTolerance.inRange(motherMass-pepMass) &&
					offsetTolerance.inRange(nTermOffset-dbHit.getNTermOffset()) &&
					offsetTolerance.inRange(cTermOffset-dbHit.getCTermOffset()))//*/
							results.add(dbHit);
			}
			
			/* Don't get it from child nodes. Each non-leaf node also has PeptideDBHit for it's position
			for (int i=0; i<AminoAcid.getIndexSize(); i++)
				if (childs[i]!=null) childs[i].getHitResult(motherMass, nTermOffset, cTermOffset, motherMassTolerance, offsetTolerance, results);
			*/
		}
	}

	public PeptideDBForPTMSearchi(){}
	
	public PeptideDBForPTMSearchi( int minQuerySeqLength, int maxQuerySeqLength, double minMass, double maxMass, int miscleavage )
	{
		this.minQuerySeqLength	= minQuerySeqLength;
		this.maxQuerySeqLength	= maxQuerySeqLength;
		this.minMass 			= minMass;
		this.maxMass			= maxMass;
		this.miscleavage		= miscleavage;
	}

	public void setQuerySeqLength(int min, int max)
	{
		this.minQuerySeqLength = min;
		this.maxQuerySeqLength = max;
	}
	
	public int construct( Iterable<Peptide> pepContainer, int minQuerySeqLength, int maxQuerySeqLength, double minMass, double maxMass, int miscleavage ) {
		this.clear();
		this.minQuerySeqLength	= minQuerySeqLength;
		this.maxQuerySeqLength	= maxQuerySeqLength;
		this.minMass 			= minMass;
		this.maxMass			= maxMass;
		this.miscleavage		= miscleavage;

        for (Peptide newPep : pepContainer)
		{
			this.put(newPep, newPep);
		}
		
		constructSuffixTree();
		
		return 0;
		
	}

	public ArrayList<PeptideDBHit> search( double motherMass, Sequence seq, double nTermOffset, double cTermOffset, DoublePair motherMassTolerance, DoublePair offsetTolerance )
	{
		ArrayList<PeptideDBHit> results = new ArrayList<>();
		
		// Wrong query
		if (seq.size()<minQuerySeqLength || seq.size()>maxQuerySeqLength) return null;
		
		// get terminals
		ArrayList<TreeNode> terminals = new ArrayList<>();
		getTerminals(seq, 0, root, terminals);
		
		for (TreeNode node : terminals)
		{
			node.xxx_getHitResult(motherMass, nTermOffset, cTermOffset, motherMassTolerance, offsetTolerance, results);			
		}
		
		return results;
	}
	
	// Get terminals of sequence
	@SuppressWarnings("DataFlowIssue")
    private void getTerminals(Sequence seq, int pos, TreeNode node, ArrayList<TreeNode> terminals )
	{
		if (node==null) return;
		if (pos==seq.size())
		{
			terminals.add(node);
			return;
		}
		if (seq.get(pos)==AminoAcid.getAminoAcid('@'))
		{
			getTerminals(seq, pos+1, node.childs[AminoAcid.getAminoAcid('I').getIndex()], terminals);
			getTerminals(seq, pos+1, node.childs[AminoAcid.getAminoAcid('L').getIndex()], terminals);			
		}
		else if (seq.get(pos)==AminoAcid.getAminoAcid('~'))	// deamidated N
		{
			getTerminals(seq, pos+1, node.childs[AminoAcid.getAminoAcid('N').getIndex()], terminals);
		}
		else if (seq.get(pos)==AminoAcid.getAminoAcid('`'))	// deamidated Q
		{
			getTerminals(seq, pos+1, node.childs[AminoAcid.getAminoAcid('Q').getIndex()], terminals);
		}
		else  
			getTerminals(seq, pos+1, node.childs[seq.get(pos).getIndex()], terminals);
	}

	
	// If there is a peptide which has same sequence as input,
	// add SrcProteinInfo to the peptide.
	// else, make new Peptide
	protected void addPeptide( Sequence seq, Protein srcProtein, int start, int end )
	{
		// Ignore sequence with to short or to long length
		if (seq.size() < minPeptideSeqLength || seq.size() > maxPeptideSeqLength) return;
		
		// Ignore sequence with too large mass
		// calculating mass every time. Might be optimized later
		if (seq.getMonoMass() > Constants.maxPeptideMass || seq.getMonoMass() < Constants.minPeptideMass) return;
	
		Peptide pep = this.get(seq);
		
		if (pep==null)
		{
			// no peptide with input sequence found, so make new one
			pep = new Peptide(srcProtein.getID(), start, end);
			pep.addAll(seq);
			this.put(pep, pep);
		}
		else
		{
			// peptide with same sequence already exists. Add SrcProteinInfo to the peptide
			pep.getSrcProteinInfo().add( new SrcProteinInfo(srcProtein.getID(), start, end) );
		}
	}
	
	protected int constructSuffixTree() {
		for (Peptide pep : this.values())
		{
			// Find all peptide sub-sequences and insert them to bucket
			double frontMass	= 0.0;
			double rearMass		= 0.0;
			for (AminoAcid aa : pep)
				rearMass += aa.getMass();
			
			for (int j=0; j<pep.size(); j++)
			{
				double		seqMass = 0.0;
				// Find a node
				TreeNode	now = root;
				int			k;
				for (k=j; k<j+maxQuerySeqLength && k<pep.size(); k++)
				{
					TreeNode next = now.childs[pep.get(k).getIndex()];
					if (next == null)
					{
						next = new TreeNode();
						now.childs[pep.get(k).getIndex()] = next;
					}
					now = next;
					
					// add mass
					seqMass += pep.get(k).getMass();

					// Insert DB Hit data to the node
					if (minQuerySeqLength<=k-j+1)
						now.bucket.add( new PeptideDBHit(pep, j, k, frontMass, rearMass-seqMass ) );
				}
				
				// calc mass
				frontMass += pep.get(j).getMass();
				rearMass -= pep.get(j).getMass();
			}
		}
		
		return 0;		
	}

	public void writeToFile(String fileName) throws IOException
	{
		ObjectOutputStream out = new ObjectOutputStream( new FileOutputStream(fileName) );
		out.writeObject(this);
		out.close();
	}
	
	public static PeptideDBForPTMSearchi readFromFile(String fileName) throws IOException, ClassNotFoundException
	{
		ObjectInputStream in = new ObjectInputStream( new FileInputStream(fileName) );
		PeptideDBForPTMSearchi pepDB = (PeptideDBForPTMSearchi)in.readObject();
		in.close();
		
		return pepDB;
	}
}