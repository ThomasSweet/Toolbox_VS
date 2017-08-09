package de.sb.toolbox.gui.swing;

import java.beans.BeanInfo;
import java.beans.IntrospectionException;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import javax.swing.table.AbstractTableModel;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.Reflections;


/**
 * Swing table model based on a list of JavaBeans.
 */
@Copyright(year = 2012, holders = "Sascha Baumeister")
public class BeanTableModel<T> extends AbstractTableModel {
	static private final long serialVersionUID = 1L;

	private final Class<T> beanClass;
	private final PropertyDescriptor[] columnDescriptors;
	private final boolean[] columnEditableFlags;
	private final List<T> beans;


	/**
	 * Creates a new instance based on introspected properties of the given bean class. Note that any indexed and write-only
	 * properties are not included within the resulting table model, as these are unsuitable for table rendering. Also note that
	 * the resulting columns are not editable by default.
	 * @param beanClass the bean class
	 * @throws NullPointerException if the given bean class is {@code null}
	 * @throws IllegalArgumentException if the given bean class cannot be introspected
	 */
	public BeanTableModel (final Class<T> beanClass) {
		this.beanClass = beanClass;
		this.columnDescriptors = introspect(beanClass);
		this.columnEditableFlags = new boolean[this.columnDescriptors.length];
		this.beans = new ArrayList<T>();
	}


	/**
	 * Returns the class of this model's row objects.
	 * @return the bean class
	 */
	public Class<T> getBeanClass () {
		return this.beanClass;
	}


	/**
	 * Returns the class of the property associated with the given column. Note that primitive types are returned as their
	 * respective wrapper classes.
	 * @param columnIndex the column index
	 * @return the column class
	 * @throws IndexOutOfBoundsException if the given column index is invalid
	 */
	@Override
	public Class<?> getColumnClass (final int columnIndex) {
		final Class<?> type = this.columnDescriptors[columnIndex].getPropertyType();
		return Reflections.boxType(type);
	}


	/**
	 * Returns the column name associated with the given index. Note that in opposition to basic table models, column names are
	 * unique for this model.
	 * @param columnIndex the column index
	 * @return the column name
	 * @throws IndexOutOfBoundsException if the given column index is invalid
	 * @see #findColumn(String)
	 */
	@Override
	public final String getColumnName (final int columnIndex) {
		return this.columnDescriptors[columnIndex].getName();
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public final int getColumnCount () {
		return this.columnDescriptors.length;
	}


	/**
	 * {@inheritDoc}
	 */
	@Override
	public final int getRowCount () {
		return this.beans.size();
	}


	/**
	 * {@inheritDoc}
	 * @throws IndexOutOfBoundsException if the given row or column index is invalid
	 */
	@Override
	public final Object getValueAt (final int rowIndex, final int columnIndex) {
		final PropertyDescriptor propertyDescriptor = this.columnDescriptors[columnIndex];
		final T bean;
		synchronized (this.beans) {
			bean = this.beans.get(rowIndex);
		}

		try {
			return propertyDescriptor.getReadMethod().invoke(bean);
		} catch (final Exception exception) {
			return null;
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IndexOutOfBoundsException if the given row or column index is invalid
	 */
	@Override
	public final void setValueAt (final Object value, final int rowIndex, final int columnIndex) {
		if (!this.isCellEditable(rowIndex, columnIndex)) return;

		final PropertyDescriptor propertyDescriptor = this.columnDescriptors[columnIndex];
		final T bean;
		synchronized (this.beans) {
			bean = this.beans.get(rowIndex);
		}

		try {
			propertyDescriptor.getWriteMethod().invoke(bean, value);
		} catch (final Exception exception) {
			return;
		}
		this.fireTableCellUpdated(rowIndex, rowIndex);
	}


	/**
	 * Returns {@code true} if the given column is editable, {@code false} otherwise. Note that columns based on read-only
	 * properties are never editable.
	 * @throws IndexOutOfBoundsException if the given row or column index is invalid
	 */
	@Override
	public boolean isCellEditable (final int rowIndex, final int columnIndex) {
		return this.isColumnEditable(columnIndex);
	}


	/**
	 * Returns {@code true} if the given column is editable, {@code false} otherwise. Note that columns based on read-only
	 * properties are never editable.
	 * @param columnIndex the column index
	 * @return {@code true} if the associated column is editable, {@code false} otherwise
	 * @throws IndexOutOfBoundsException if the given column index is invalid
	 */
	public boolean isColumnEditable (final int columnIndex) {
		return this.columnEditableFlags[columnIndex];
	}


	/**
	 * Defines if the given column is editable or not. Note that columns based on read-only properties can never be editable,
	 * therefore attempts to set them editable are silently ignored.
	 * @param columnIndex the column index
	 * @param editable {@code true} if the associated column is editable, {@code false} otherwise
	 * @throws IndexOutOfBoundsException if the given column index is invalid
	 */
	public void setColumnEditable (final int columnIndex, final boolean editable) {
		if (this.columnDescriptors[columnIndex].getWriteMethod() == null) return;

		this.columnEditableFlags[columnIndex] = editable;
	}


	/**
	 * Adds the given beans to this model's rows, and fires a change event.
	 * @param beans the row beans to be added
	 * @throws NullPointerException if the given bean array or any of it's elements is {@code null}
	 * @throws IllegalArgumentException if the model's bean class cannot be assigned from one of the given bean's classes
	 */
	public void addRows (final T[] beans) {
		if (beans.length == 0) return;

		final int rowIndex;
		synchronized (this.beans) {
			rowIndex = this.beans.size();

			for (final T bean : beans) {
				if (!this.beanClass.isAssignableFrom(bean.getClass())) throw new IllegalArgumentException();

				this.beans.add(bean);
			}
		}
		this.fireTableRowsInserted(rowIndex, rowIndex + beans.length - 1);
	}


	/**
	 * Adds a bean to this model's rows, and fires a change event.
	 * @param bean the row bean to be added
	 * @throws NullPointerException if the given bean is {@code null}
	 * @throws IllegalArgumentException if the model's bean class cannot be assigned from the given bean's class
	 */
	public void addRow (final T bean) {
		if (!this.beanClass.isAssignableFrom(bean.getClass())) throw new IllegalArgumentException();

		final int rowIndex;
		synchronized (this.beans) {
			rowIndex = this.beans.size();
			this.beans.add(bean);
		}
		this.fireTableRowsInserted(rowIndex, rowIndex);
	}


	/**
	 * Inserts a bean into this model's rows at the given row index, and fires a change event.
	 * @param rowIndex the row index
	 * @param bean the row bean to be added
	 * @throws NullPointerException if the given bean is {@code null}
	 * @throws IllegalArgumentException if the model's bean class cannot be assigned from the given bean's class
	 * @throws ArrayIndexOutOfBoundsException if the given row index is invalid
	 */
	public void insertRow (final int rowIndex, final T bean) {
		if (!this.beanClass.isAssignableFrom(bean.getClass())) throw new IllegalArgumentException();

		synchronized (this.beans) {
			this.beans.add(rowIndex, bean);
		}

		this.fireTableRowsInserted(rowIndex, rowIndex);
	}


	/**
	 * Returns a row bean from this model's rows. Note that the beans should never be modified directly, as this does not fire
	 * the necessary change events to reflect the updates within the table!
	 * @return the row bean at the given row index
	 * @throws ArrayIndexOutOfBoundsException if the given row index is invalid
	 */
	public T getRow (final int rowIndex) {
		synchronized (this.beans) {
			return this.beans.get(rowIndex);
		}
	}


	/**
	 * Returns this model's beans as an array. Note that the beans should never be modified directly, as this does not fire the
	 * necessary change events to reflect the updates within the table!
	 * @return the model's row beans
	 */
	public T[] getRows () {
		@SuppressWarnings("unchecked")
		final T[] template = (T[]) Array.newInstance(this.beanClass, 0);
		return this.beans.toArray(template);
	}


	/**
	 * Moves the inclusive row bean range defined by {@code fromIndex} and {@code toIndex} up or down this model, and fires a
	 * change event. If the given distance is negative, the range is moved upwards by {@code -distance} rows. Otherwise, the
	 * range is moved downwards by {@code distance} rows.
	 * @param fromIndex the from index
	 * @param toIndex the to index
	 * @param distance the distance and direction the range shall be moved
	 * @throws IndexOutOfBoundsException if one of the given indices is out of range
	 */
	public void moveRows (final int fromIndex, final int toIndex, final int distance) {
		if (distance == 0) return;

		synchronized (this.beans) {
			final List<T> beanRange = this.beans.subList(fromIndex, toIndex);
			final List<T> elements = new ArrayList<>(beanRange);
			beanRange.clear();
			this.beans.addAll(fromIndex + distance, elements);
		}

		final int firstRow = distance < 0 ? fromIndex + distance : fromIndex;
		final int lastRow = distance < 0 ? toIndex : toIndex + distance;
		this.fireTableRowsUpdated(firstRow, lastRow);
	}


	/**
	 * Removes a row bean from this model at the given row index, and fires a change event.
	 * @param rowIndex the row index
	 * @throws ArrayIndexOutOfBoundsException if the given row index is invalid
	 */
	public void removeRow (final int rowIndex) {
		synchronized (this.beans) {
			this.beans.remove(rowIndex);
		}

		this.fireTableRowsDeleted(rowIndex, rowIndex);
	}


	/**
	 * Removes all row beans from this model, and fires a change event.
	 */
	public void removeRows () {
		final int rowCount;
		synchronized (this.beans) {
			rowCount = this.beans.size();
			this.beans.clear();
		}

		if (rowCount > 0) {
			this.fireTableRowsDeleted(0, rowCount - 1);
		}
	}


	/**
	 * Introspects the given bean class, and returns descriptors for all properties that can be accessed as part of a table
	 * model. This excludes any indexed and write-only properties.
	 * @param beanClass the bean class
	 * @return the property descriptors
	 */
	static private PropertyDescriptor[] introspect (final Class<?> beanClass) {
		final BeanInfo beanInfo;
		try {
			beanInfo = Introspector.getBeanInfo(beanClass);
		} catch (final IntrospectionException exception) {
			throw new IllegalArgumentException(exception);
		}

		final List<PropertyDescriptor> result = new ArrayList<>();
		for (final PropertyDescriptor descriptor : beanInfo.getPropertyDescriptors()) {
			if (descriptor.getPropertyType() != null & descriptor.getReadMethod() != null) {
				result.add(descriptor);
			}
		}
		return result.toArray(new PropertyDescriptor[0]);
	}
}