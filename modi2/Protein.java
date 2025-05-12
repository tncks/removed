package modi;

@SuppressWarnings("unused")
public class Protein extends Sequence {
	final String 		annotation;
	String 		name;
	final int			id;		// index in protein container
	
	public String			getAnnotaion()	{ return annotation; }
	public String			getName()		{ return name; }
	public int				getID()			{ return id; }
	
	public Protein( String h, int id ){
		int i, cut= 0;
		for( i=0; i<h.length(); i++ ){
			if( h.charAt(i) == '|' || h.charAt(i) == ':' || h.charAt(i) == '_' ) cut++;			
			if( cut==1 || h.charAt(i) == ' ') break;
		}				
		annotation = h.substring(1, i);
		this.id				= id;
	}
	
	public String toString()
	{
		return "PROTEIN("+annotation+")="+super.toString();
	}
	
	public String getSequence()
	{
		return super.toString();
	}
	
	public Protein getReversed()
	{
		Protein reverse = new Protein(annotation+"r", id);
		for(int i=0; i<this.size(); i++)
		{
			reverse.add(this.get(this.size()-1-i));
		}
		return reverse;
	}
// test by timebird
}
