package structures;

import java.util.Arrays;
import java.util.Random;

import shared.Tools;
import stream.KillSwitch;

/**
 * Like IntHashMap, but uses a power-of-2 length to avoid modulo operations.
 * @author Brian Bushnell
 * @date June 8, 2017
 *
 */
public final class IntHashMapBinary extends AbstractIntHashMap{	
	
	public static void main(String[] args){
		IntHashMapBinary set=new IntHashMapBinary(32, 0.7f);
		test(set);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public IntHashMapBinary(){
		this(256);
	}
	
	public IntHashMapBinary(int initialSize){
		this(initialSize, 0.7f);
	}
	
	public IntHashMapBinary(int initialSize, float loadFactor_){
		if(Integer.bitCount(initialSize)>1){
			int zeros=Integer.numberOfLeadingZeros(initialSize);
			if(zeros<2){zeros=2;}
			initialSize=1<<(32-zeros);
		}
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
		
		final int limit=koys.length, initial=key&modulus;
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
		
		final int limit=koys.length, initial=key&modulus;
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
		resize(Tools.max(2, modulus+1)*2L);
	}
	
	private final void resize(final long size2){
		assert(size2>size) : size+", "+size2;
		assert(Long.bitCount(size2)==1) : size+", "+size2;
		long newModulus=size2-1;
		assert(newModulus>0) : newModulus;
		assert(newModulus+extra<Integer.MAX_VALUE) : "Overflow";
		assert(newModulus>modulus);
		
//		System.err.println("size="+size+", modulus="+modulus+" -> size2="+size2+", modulus2="+newModulus);
		
		modulus=(int)newModulus;
		
		final int size3=(int)(newModulus+extra);
		sizeLimit=(int)(newModulus*loadFactor);
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
//	private static final int MASK=Integer.MAX_VALUE;
	private static final int MINMASK=Integer.MIN_VALUE;
	
	private static final int extra=10;
	
}
