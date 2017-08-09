package de.sb.toolbox.util;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.util.AbstractSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.Set;
import de.sb.toolbox.Copyright;


/**
 * This is a wrapper class enriching underlying sets with event callback functionality. Instances notify their registered event
 * handlers of every change event within the underlying set (via {@link Set#add(Object)}, {@link Set#remove(Object)}, or
 * {@link Set#clear()} messages). This includes modifications caused by associated iterators. Note that event notification
 * doesn't happen if the underlying set (the one provided with the constructor) is modified directly.
 * @param <E> the element type
 */
@Copyright(year = 2010, holders = "Sascha Baumeister")
public class CallbackSet<E> extends AbstractSet<E>implements Set<E> {

	private final Set<E> delegateSet;
	private final Set<VetoableChangeListener> listeners;


	/**
	 * Creates a new instance based on an empty hash map.
	 */
	public CallbackSet () {
		this(new HashSet<E>());
	}


	/**
	 * Creates a new instance based on the delegate set. Note that no put events are spawned for the elements already existing
	 * within the delegate set.
	 * @param delegateSet the underlying set
	 * @see Set#addAll(Collection)
	 */
	public CallbackSet (final Set<E> delegateSet) {
		this.delegateSet = delegateSet;
		this.listeners = Collections.synchronizedSet(new HashSet<VetoableChangeListener>());
	}


	/**
	 * Returns the registered event listeners.
	 * @return the event listeners
	 */
	public Set<VetoableChangeListener> getListeners () {
		return this.listeners;
	}


	/**
	 * Returns the underlying delegate set. Note that event notifications don't take place when modifying the delegate set
	 * directly, which is the intended use-case for this method.
	 * @return the delegate set
	 */
	public Set<E> getDelegateSet () {
		return this.delegateSet;
	}


	/**
	 * {@inheritDoc}
	 * @throws ClassCastException {@inheritDoc}
	 * @throws NullPointerException {@inheritDoc}
	 */
	@Override
	public boolean contains (final Object value) {
		return this.delegateSet.contains(value);
	}


	/**
	 * {@inheritDoc}
	 * @throws UnsupportedOperationException {@inheritDoc}
	 * @throws NullPointerException {@inheritDoc}
	 * @throws ClassCastException {@inheritDoc}
	 * @throws IllegalArgumentException {@inheritDoc}
	 * @throws IllegalStateException if a listener vetoes the change
	 */
	@Override
	public boolean add (final E value) {
		if (this.delegateSet.contains(value)) return false;

		this.fireChangeEvent("add", null, value);
		return this.delegateSet.add(value);
	}


	/**
	 * {@inheritDoc}
	 * @throws UnsupportedOperationException {@inheritDoc}
	 * @throws NullPointerException {@inheritDoc}
	 * @throws ClassCastException {@inheritDoc}
	 * @throws IllegalStateException if a listener vetoes the change
	 */
	@Override
	public boolean remove (final Object value) {
		try {
			final E element = this.find(value);
			this.fireChangeEvent("remove", element, null);
			return this.delegateSet.remove(value);
		} catch (final NoSuchElementException exception) {
			return false;
		}
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public Iterator<E> iterator () {
		return new CallbackIterator();
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public int size () {
		return this.delegateSet.size();
	}


	/**
	 * Returns the set element that is equal to the given value.
	 * @param value the value
	 * @return the element
	 * @throws NoSuchElementException if the given value is not contained within this set
	 */
	private E find (final Object value) {
		if (this.delegateSet.contains(value)) {
			for (final E element : this.delegateSet) {
				if (Objects.equals(element, value)) return element;
			}
		}

		throw new NoSuchElementException();
	}


	/**
	 * Notifies all listeners of an impeding change. If one of the listeners throws a {@link PropertyVetoException} during
	 * notification, this set will not perform the change.
	 * @param key the key
	 * @param oldValue the old value
	 * @param newValue the new value
	 * @see VetoableChangeListener#vetoableChange(PropertyChangeEvent)
	 * @throws IllegalStateException if one of the listeners vetoes the change
	 */
	private void fireChangeEvent (final String key, final E oldValue, final E newValue) {
		final PropertyChangeEvent event = new PropertyChangeEvent(this, key, oldValue, newValue);
		final VetoableChangeListener[] listeners = CallbackSet.this.listeners.toArray(new VetoableChangeListener[0]);
		for (final VetoableChangeListener listener : listeners) {
			try {
				listener.vetoableChange(event);
			} catch (final PropertyVetoException exception) {
				throw new IllegalStateException(exception);
			}
		}
	}



	/**
	 * Inner class defining a callback iterator.
	 */
	private final class CallbackIterator implements Iterator<E> {
		private final Iterator<E> iterator = CallbackSet.this.delegateSet.iterator();
		private E currentElement = null;


		@Override
		public boolean hasNext () {
			return this.iterator.hasNext();
		}


		@Override
		public E next () {
			return this.currentElement = this.iterator.next();
		}


		@Override
		public void remove () {
			CallbackSet.this.fireChangeEvent("remove", this.currentElement, null);
			this.iterator.remove();
		}
	}
}