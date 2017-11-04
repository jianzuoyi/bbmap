package structures;

import java.util.Arrays;
import java.util.Random;

import shared.Primes;
import shared.Tools;
import stream.KillSwitch;

/**
 * @author Brian Bushnell
 * @date June 7, 2017
 *
 */
public final class IntHashMap extends AbstractIntHashMap{	
	
	public static void main(String[] args){
		IntHashMap set=new IntHashMap(20, 0.7f);
		test(set);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public IntHashMap(){
		this(256);
	}
	
	public IntHashMap(int initialSize){
		this(initialSize, 0.7f);
	}
	
	public IntHashMap(int initialSize, float loadFactor_){
		invalid=randy.nextInt()|MINMASK;
		assert(invalid<0);
		assert(initialSize>0);
		assert(loadFactor_>0 && loadFactor_<1);
		loadFactor=Tools.mid(0.25f, loadFactor_, 0.90f);
		resize(initialSize);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Public Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	public void clear(){
		if(size<1){return;}
		Arrays.fill(koys, invalid);
		Arrays.fill(volues, 0);
		size=0;
	}
	
	@Override
	public int get(int key){
		int cell=findCell(key);
		return cell<0 ? -1 : volues[cell];
	}
	
	@Override
	public int put(int key, int value){return set(key, value);}
	
	@Override
	public int set(int key, int value){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=volues[cell];
		volues[cell]=value;
		if(koys[cell]==invalid){
			koys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return oldV;
	}
	
	@Override
	public int increment(int key){
		return increment(key, 1);
	}
	
	@Override
	public int increment(int key, int incr){
		if(key==invalid){resetInvalid();}
		final int cell=findCellOrEmpty(key);
		final int oldV=volues[cell];
		final int value=oldV+incr;
		volues[cell]=value;
		volues[cell]=Tools.min(Integer.MAX_VALUE, value);
		if(koys[cell]==invalid){
			koys[cell]=key;
			size++;
			if(size>sizeLimit){resize();}
		}
//		assert(get(key)==value);//TODO: slow
		return value;
	}
	
	@Override
	public boolean remove(int key){
		if(key==invalid){return false;}
		final int cell=findCell(key);
		if(cell<0){return false;}
		assert(koys[cell]==key);
		koys[cell]=invalid;
		volues[cell]=0;
		size--;
		
		rehashFrom(cell);
		return true;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Private Methods       ----------------*/
	/*--------------------------------------------------------------*/
	
	private void rehashFrom(int initial){
		if(size<1){return;}
		final int limit=koys.length;
		for(int cell=initial+1; cell<limit; cell++){
			final int key=koys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
		for(int cell=0; cell<initial; cell++){
			final int key=koys[cell];
			if(key==invalid){return;}
			rehashCell(cell);
		}
	}
	
	private boolean rehashCell(final int cell){
		final int key=koys[cell];
		final int value=volues[cell];
		assert(key!=invalid);
		if(key==invalid){resetInvalid();}
		final int dest=findCellOrEmpty(key);
		if(cell==dest){return false;}
		assert(koys[dest]==invalid);
		koys[cell]=invalid;
		volues[cell]=0;
		koys[dest]=key;
		volues[dest]=value;
		
		return true;
	}
	
	private void resetInvalid(){
		final int old=invalid;
		int x=invalid;
		while(x==old || contains(x)){x=randy.nextInt()|MINMASK;}
		assert(x<0);
		invalid=x;
		for(int i=0; i<koys.length; i++){
			if(koys[i]==old){
				koys[i]=invalid;
//				assert(volues[i]==0); //TODO: slow
			}
		}
	}

	@Override
	int findCell(final int key){
		if(key==invalid){return -1;}
		
		final int limit=koys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final int x=koys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=koys[cell];
			if(x==key){return cell;}
			if(x==invalid){return -1;}
		}
		return -1;
	}
	
	private int findCellOrEmpty(final int key){
		assert(key!=invalid) : "Collision - this should have been intercepted.";
		
		final int limit=koys.length, initial=(int)((key&MASK)%modulus);
		for(int cell=initial; cell<limit; cell++){
			final int x=koys[cell];
			if(x==key || x==invalid){return cell;}
		}
		for(int cell=0; cell<initial; cell++){
			final int x=koys[cell];
			if(x==key || x==invalid){return cell;}
		}
		throw new RuntimeException("No empty cells - size="+size+", limit="+limit);
	}
	
	private final void resize(){
		assert(size>=sizeLimit);
		resize(koys.length*2L+1);
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		long newPrime=Primes.primeAtLeast(size2);
		if(newPrime+extra>Integer.MAX_VALUE){
			newPrime=Primes.primeAtMost(Integer.MAX_VALUE-extra);
		}
		assert(newPrime>modulus) : "Overflow: "+size+", "+size2+", "+modulus+", "+newPrime;
		modulus=(int)newPrime;
		
		final int size3=(int)(newPrime+extra);
		sizeLimit=(int)(modulus*loadFactor);
		final int[] oldK=koys;
		final int[] oldV=volues;
		koys=KillSwitch.allocInt1D(size3);
		volues=KillSwitch.allocInt1D(size3);
		Arrays.fill(koys, invalid);
		
//		System.err.println("Resizing "+(old==null ? "null" : ""+old.length)+" to "+size3);
		
		if(size<1){return;}
		
		size=0;
		for(int i=0; i<oldK.length; i++){
			final int k=oldK[i], v=oldV[i];
			if(k!=invalid){
				set(k, v);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Getters           ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	public int[] keys(){return koys;}
	
	@Override
	public int[] values(){return volues;}
	
	@Override
	public int invalid(){return invalid;}
	
	@Override
	public int size(){return size;}
	
	@Override
	public boolean isEmpty(){return size==0;}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int[] koys;
	private int[] volues;
	private int size=0;
	/** Value for empty cells */
	private int invalid;
	private int modulus;
	private int sizeLimit;
	private final float loadFactor;
	
	private static final Random randy=new Random(1);
	private static final int MASK=Integer.MAX_VALUE;
	private static final int MINMASK=Integer.MIN_VALUE;
	
	private static final int extra=10;
	
}
