package ClassicalStat;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

//Similar a <see cref="System.Collections.Generic.SortedDictionary{TKey, TValue}"/> pero con acceso subindicado en O(1). Implementa <see cref="IDictionary{TKey, TValue}"/>

public class ECCSortedDictionary<TKey, TValue> implements Map<TKey, TValue> {

	private HashMap<TKey, TValue> internaldict = null;
	private TreeMap<TKey, Object> sortedkeys = null;

	/** 
	 Constructor.
	*/
	public ECCSortedDictionary() {
		internaldict = new HashMap<TKey, TValue>();
		sortedkeys = new TreeMap<TKey, Object>();
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.Add"/>
	*/
	public final void Add(TKey key, TValue value)	{
		internaldict.put(key, value);
		sortedkeys.put(key, null);
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.ContainsKey"/>
	*/
	public final boolean containsKey(Object objectKey)	{
		TKey key = (TKey)objectKey;
		return internaldict.containsKey(key);
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.Keys"/>
	*/
	public final Object[] getKeys() {
		return sortedkeys.keySet().toArray();
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.Remove"/>
	*/
	public final void Remove(TKey key)	{
		sortedkeys.remove(key); 
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.TryGetValue"/>
	*/
	public final boolean TryGetValue(TKey key, DotNetHelpers.RefObject<TValue> value) {
		return (internaldict.containsKey(key) ? (value.argValue = internaldict.get(key)) == value.argValue : false);
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.Values"/>
	*/
	public final Collection<TValue> values() {
		return internaldict.values();
	}

	/** 
	 Implementa <see cref="IDictionary{TKey, TValue}.Item"/>
	*/
	public final TValue get(Object objectKey)	{
		TKey key = (TKey)objectKey;
		return internaldict.get(key);
	}

	public final void setItem(TKey key, TValue value)	{
		internaldict.put(key, value);
		sortedkeys.put(key, value);
	}

	/** 
	 *NO* Implementa <see cref="IDictionary{TKey, TValue}.Add"/> de <see cref="KeyValuePair{TKey, TValue}"/>
	*/

	public final void Add(Map.Entry<TKey, TValue> item) {
		throw new UnsupportedOperationException();
	}

	/** 
	 Implementa <see cref="ICollection{T}.Clear()"/>
	*/
	public final void clear()	{
		internaldict.clear();
		sortedkeys.clear();
	}

	/** 
	 *NO* Implementa <see cref="ICollection{T}.Contains"/> de <see cref="KeyValuePair{TKey, TValue}"/>
	*/
	public final boolean contains(Object objectValue)	{
		Map.Entry<TKey, TValue> item = (Map.Entry<TKey, TValue>)objectValue;
		throw new UnsupportedOperationException();
	}

	/** 
	 *NO* Implementa <see cref="ICollection{T}.CopyTo"/> de <see cref="KeyValuePair{TKey, TValue}"/>
	*/
	public final void CopyTo(Map.Entry<TKey, TValue>[] array, int arrayIndex) {
		throw new UnsupportedOperationException();
	}

	/** 
	 Implementa <see cref="ICollection{T}.Count"/>
	*/
	public final int size() {
		return internaldict.size();
	}

	/** 
	 *NO* Implementa <see cref="ICollection{T}.IsReadOnly"/>
	*/
	public final boolean getIsReadOnly()	{
		throw new UnsupportedOperationException();
	}

	/** 
	 *NO* Implementa <see cref="IDictionary{TKey, TValue}.Remove"/> de <see cref="KeyValuePair{TKey, TValue}"/>
	*/
	public final boolean Remove(Map.Entry<TKey, TValue> item) {
		throw new UnsupportedOperationException();
	}

	/** 
	 *NO* Implementa <see cref="IEnumerable{T}.GetEnumerator"/> de <see cref="KeyValuePair{TKey, TValue}"/>
	*/
	public final Iterator<Map.Entry<TKey, TValue>> iterator()	{
		throw new UnsupportedOperationException();
	}

	public final Iterator GetEnumerator()	{
		return internaldict.entrySet().iterator();
	}

	@Override
	public boolean containsValue(Object value) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<java.util.Map.Entry<TKey, TValue>> entrySet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public Set<TKey> keySet() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TValue put(TKey key, TValue value) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void putAll(Map<? extends TKey, ? extends TValue> m) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public TValue remove(Object key) {
		// TODO Auto-generated method stub
		return null;
	}
}
