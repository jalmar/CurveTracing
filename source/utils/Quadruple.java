package utils;

/**
 *	Note: class implementation adapted from http://stackoverflow.com/questions/156275/what-is-the-equivalent-of-the-c-pairl-r-in-java and http://stackoverflow.com/questions/2670982/using-pairs-or-2-tuples-in-java
 */
public class Quadruple<A, B, C, D>
{
	// members
	private A first;
	private B second;
	private C third;
	private D fourth;
	
	// *************************************************************************
	
	// constructor
	public Quadruple()
	{
		super();
		this.first = null;
		this.second = null;
		this.third = null;
		this.fourth = null;
	}
	
	// NOTE: cannot do single argument constructor in case type A == type B
	
	public Quadruple(A first, B second, C third, D fourth)
	{
		super();
		this.first = first;
		this.second = second;
		this.third = third;
		this.fourth = fourth;
	}
	
	// *************************************************************************
	
	public A getFirst()
	{
		return this.first;
	}
	
	public void setFirst(A first)
	{
		this.first = first;
	}
	
	public B getSecond()
	{
		return this.second;
	}
	
	public void setSecond(B second)
	{
		this.second = second;
	}
	
	public C getThird()
	{
		return this.third;
	}
	
	public void setThird(C third)
	{
		this.third = third;
	}
	
	public D getFourth()
	{
		return this.fourth;
	}
	
	public void setFourth(D fourth)
	{
		this.fourth = fourth;
	}
	
	// *************************************************************************
	
	@Override
	public int hashCode() {
		int hashFirst = this.first != null ? this.first.hashCode() : 0;
		int hashSecond = this.second != null ? this.second.hashCode() : 0;
		int hashThird = this.third != null ? this.third.hashCode() : 0;
		int hashFourth = this.fourth != null ? this.fourth.hashCode() : 0;
		
		return (hashFirst + hashSecond + hashThird + hashFourth) * hashFourth + (hashFirst + hashSecond + hashThird) * hashThird + (hashFirst + hashSecond) * hashSecond + hashFirst;
		
		//final int prime = 31;
		//int result = 1;
		//result = prime * result + ((x == null) ? 0 : x.hashCode());
		//result = prime * result + ((y == null) ? 0 : y.hashCode());
		//return result;
	}
	
	@Override
	public boolean equals(Object other)
	{
		if (other instanceof Quadruple) {
			Quadruple otherQuadruple = (Quadruple) other;
			return 
			((  this.first == otherQuadruple.first ||
				( this.first != null && otherQuadruple.first != null &&
				  this.first.equals(otherQuadruple.first))) &&
			 (	this.second == otherQuadruple.second ||
				( this.second != null && otherQuadruple.second != null &&
				  this.second.equals(otherQuadruple.second))) &&
			 (	this.third == otherQuadruple.third ||
				( this.third != null && otherQuadruple.third != null &&
				  this.third.equals(otherQuadruple.third))) &&
			 (	this.fourth == otherQuadruple.fourth ||
				( this.fourth != null && otherQuadruple.fourth != null &&
				  this.fourth.equals(otherQuadruple.fourth)))
			);
		}
		
		return false;
	}
	
	public String toString()
	{ 
		return "(" + this.first + ", " + this.second + ", " + this.third + ", " + this.fourth + ")"; 
	}
}
