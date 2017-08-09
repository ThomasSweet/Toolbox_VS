package de.sb.toolbox.util;

import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.SQLFeatureNotSupportedException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import de.sb.toolbox.Copyright;


/**
 * Instances of subtypes of this class adapt a JDBC result set (cursors) into a {@link ListIterator} of row elements. Subclasses
 * are free to choose a suitable element type for the rows; the default implementations model rows as either maps associating
 * column labels with cell values, or as lists of cell values.
 */
@Copyright(year = 2016, holders = "Sascha Baumeister")
public abstract class ResultSetIterator<T> implements ListIterator<T> {
	private final ResultSet cursor;
	private final List<String> columnLabels;


	public ResultSetIterator (final ResultSet cursor) throws NullPointerException, IllegalStateException {
		final String[] columnLabels;
		try {
			final ResultSetMetaData metaData = cursor.getMetaData();
			columnLabels = new String[metaData.getColumnCount()];
			for (int index = 0; index < columnLabels.length; ++index) {
				columnLabels[index] = metaData.getColumnLabel(index + 1);
			}
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}

		this.cursor = cursor;
		this.columnLabels = Collections.unmodifiableList(Arrays.asList(columnLabels));
	}


	/**
	 * Returns the table cursor.
	 * @return the cursor
	 */
	public ResultSet getCursor () {
		return this.cursor;
	}


	/**
	 * Returns the table's column labels.
	 * @return the column labels as an unmodifiable list
	 */
	public List<String> getColumnLabels () {
		return this.columnLabels;
	}


	/**
	 * Returns the current table row.
	 * @return the table row
	 * @throws SQLException if a database access error occurs or this method is called on a closed result set
	 */
	protected abstract T getRow () throws SQLException;


	/**
	 * {@inheritDoc}
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument is illegal
	 * @throws SQLFeatureNotSupportedException if the JDBC driver does not support this operation
	 * @throws SQLException if there is a database related problem
	 */
	protected abstract void setRow (T row) throws NullPointerException, IllegalArgumentException, SQLException;


	/**
	 * Returns the last element returned by {@link #next} or {@link #previous}.
	 * @throws IllegalStateException if there is a database related problem
	 */
	public final T get () throws IllegalStateException {
		try {
			return this.getRow();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument is illegal
	 * @throws IllegalStateException {@inheritDoc}, or if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code set} operation is not supported by this list iterator
	 */
	public final void set (final T row) throws NullPointerException, IllegalArgumentException, IllegalStateException, UnsupportedOperationException {
		try {
			this.setRow(row);
			this.cursor.updateRow();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code hasPrevious} operation is not supported by this list iterator
	 */
	public final boolean hasPrevious () throws IllegalStateException {
		try {
			return !this.cursor.isBeforeFirst() && !this.cursor.isFirst();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code hasNext} operation is not supported by this list iterator
	 */
	public final boolean hasNext () throws IllegalStateException, UnsupportedOperationException {
		try {
			return !this.cursor.isLast() && !this.cursor.isAfterLast();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws NoSuchElementException {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code previous} operation is not supported by this list iterator
	 */
	public final T previous () throws NoSuchElementException, IllegalStateException, UnsupportedOperationException {
		try {
			if (this.cursor.previous()) return this.getRow();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}

		throw new NoSuchElementException();
	}


	/**
	 * {@inheritDoc}
	 * @throws NoSuchElementException {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 */
	public final T next () throws NoSuchElementException, IllegalStateException {
		try {
			if (this.cursor.next()) return this.getRow();
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}

		throw new NoSuchElementException();
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code previousIndex} operation is not supported by this list iterator
	 */
	public final int previousIndex () throws IllegalStateException, UnsupportedOperationException {
		try {
			return this.cursor.isBeforeFirst() ? -1 : this.cursor.getRow() - 2;
		} catch (final SQLFeatureNotSupportedException exception) {
			throw new UnsupportedOperationException(exception);
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException if the {@code nextIndex} operation is not supported by this list iterator
	 */
	public final int nextIndex () throws IllegalStateException, UnsupportedOperationException {
		try {
			return this.cursor.isAfterLast() ? this.cursor.getRow() - 1 : this.cursor.getRow();
		} catch (final SQLFeatureNotSupportedException exception) {
			throw new UnsupportedOperationException(exception);
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException {@inheritDoc}, or if there is a database related problem
	 * @throws UnsupportedOperationException {@inheritDoc}
	 */
	public final void remove () {
		try {
			this.cursor.deleteRow();
		} catch (final SQLFeatureNotSupportedException exception) {
			throw new UnsupportedOperationException(exception);
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if there is a database related problem
	 * @throws UnsupportedOperationException {@inheritDoc}
	 */
	public final void add (final T row) throws IllegalStateException, UnsupportedOperationException {
		try {
			this.cursor.moveToInsertRow();
			try {
				this.setRow(row);
				this.cursor.insertRow();
			} finally {
				this.cursor.moveToCurrentRow();
			}
		} catch (final SQLFeatureNotSupportedException exception) {
			throw new UnsupportedOperationException(exception);
		} catch (final SQLException exception) {
			throw new IllegalStateException(exception);
		}
	}


	/**
	 * Returns a new result set iterator that exposes rows as maps associating column labels with cell values.
	 * @param cursor the cursor
	 * @param sortedColumns whether the columns shall be sorted or not
	 * @return the map based result set iterator
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalArgumentException if the given cursor's column labels are not unique
	 * @throws IllegalStateException if a database access error occurs or the given cursor is closed
	 */
	static public ResultSetIterator<Map<String,Object>> newInstance (final ResultSet cursor, final boolean sortedColumns) throws NullPointerException, IllegalArgumentException, IllegalStateException {
		final ResultSetIterator<Map<String,Object>> iterator = new ResultSetIterator<Map<String,Object>>(cursor) {
			protected Map<String,Object> getRow () throws SQLException {
				final Map<String,Object> row = sortedColumns ? new TreeMap<>() : new HashMap<>();
				for (final String columnLabel : this.getColumnLabels()) {
					row.put(columnLabel, this.getCursor().getObject(columnLabel));
				}
				return row;
			}


			protected void setRow (final Map<String,Object> row) throws NullPointerException, SQLFeatureNotSupportedException, SQLException {
				for (final Map.Entry<String,?> cell : row.entrySet()) {
					this.getCursor().updateObject(cell.getKey(), cell.getValue());
				}
			}
		};

		if (iterator.getColumnLabels().size() == new HashSet<>(iterator.getColumnLabels()).size()) return iterator;
		throw new IllegalArgumentException();
	}


	/**
	 * Returns a new result set iterator that exposes rows as lists of cell values.
	 * @param cursor the cursor
	 * @return the list based result set iterator
	 * @throws NullPointerException if any of the given arguments is {@code null}
	 * @throws IllegalStateException if a database access error occurs or the given cursor is closed
	 */
	static public ResultSetIterator<List<Object>> newInstance (final ResultSet cursor) throws NullPointerException, IllegalStateException {
		return new ResultSetIterator<List<Object>>(cursor) {
			protected List<Object> getRow () throws SQLException {
				final List<Object> row = Arrays.asList(new Object[this.getColumnLabels().size()]);
				for (int columnIndex = 0, stop = row.size(); columnIndex < stop; ++columnIndex) {
					row.set(columnIndex, this.getCursor().getObject(columnIndex + 1));
				}
				return row;
			}


			protected void setRow (final List<Object> row) throws NullPointerException, SQLFeatureNotSupportedException, SQLException {
				for (int columnIndex = 0, stop = row.size(); columnIndex < stop; ++columnIndex) {
					this.getCursor().updateObject(columnIndex + 1, row.get(columnIndex));
				}
			}
		};
	}
}