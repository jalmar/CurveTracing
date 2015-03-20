
package utils;

/**
 *	Note: class implementation adapted from http://stackoverflow.com/questions/156275/what-is-the-equivalent-of-the-c-pairl-r-in-java and http://stackoverflow.com/questions/2670982/using-pairs-or-2-tuples-in-java
 */
public class Tuple<A, B>
{
	// members
	public A first;
	public B second;
	
	// *************************************************************************
	
	// constructor
	public Tuple()
	{
		super();
		this.first = null;
		this.second = null;
	}
	
	// NOTE: cannot do single argument constructor in case type A == type B
	
	public Tuple(A first, B second)
	{
		super();
		this.first = first;
		this.second = second;
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
	
	// *************************************************************************
	
	@Override
	public int hashCode() {
		int hashFirst = this.first != null ? this.first.hashCode() : 0;
		int hashSecond = this.second != null ? this.second.hashCode() : 0;
		return (hashFirst + hashSecond) * hashSecond + hashFirst;
		
		//final int prime = 31;
		//int result = 1;
		//result = prime * result + ((x == null) ? 0 : x.hashCode());
		//result = prime * result + ((y == null) ? 0 : y.hashCode());
		//return result;
	}
	
	@Override
	public boolean equals(Object other)
	{
		if (other instanceof Tuple) {
			Tuple otherTuple = (Tuple) other;
			return 
			((  this.first == otherTuple.first ||
				( this.first != null && otherTuple.first != null &&
				  this.first.equals(otherTuple.first))) &&
			 (	this.second == otherTuple.second ||
				( this.second != null && otherTuple.second != null &&
				  this.second.equals(otherTuple.second))) );
		}
		
		return false;
	}
	
	public String toString()
	{ 
		return "(" + this.first + ", " + this.second + ")"; 
	}
}
