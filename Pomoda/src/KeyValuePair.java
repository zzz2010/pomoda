import java.io.Serializable;
import java.util.Map.Entry;


public class KeyValuePair<K,V> implements Entry<K, V>, Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = -3225161899777880782L;
	public KeyValuePair(K key, V value) {
		super();
		this.key = key;
		this.value = value;
	}

	K key;
	V value;
	@Override
	public K getKey() {
		// TODO Auto-generated method stub
		return key;
	}

	@Override
	public V getValue() {
		// TODO Auto-generated method stub
		return value;
	}

	@Override
	public V setValue(V arg0) {
		// TODO Auto-generated method stub
		value=arg0;
		return arg0;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return key.toString()+":"+value.toString();
	}
	
	

}
