package de.sb.toolbox.gui.swing;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.text.Format;
import java.util.EventObject;
import java.util.Objects;
import javax.swing.AbstractCellEditor;
import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import javax.swing.plaf.synth.SynthContext;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import javax.swing.tree.TreeCellEditor;
import javax.swing.tree.TreeCellRenderer;
import de.sb.toolbox.Copyright;
import de.sb.toolbox.Reflections;


/**
 * Instances of this class work as combined table/tree cell renderer and editor, based on an underlying format implementation.
 * @see javax.swing.plaf.synth.SynthTableUI#paintCell(SynthContext, Graphics, Rectangle, int, int)
 */
@Copyright(year = 2012, holders = "Sascha Baumeister")
public class TableCellFormater extends AbstractCellEditor implements TableCellRenderer, TreeCellRenderer, TableCellEditor, TreeCellEditor {
	static private final long serialVersionUID = 1L;
	static private final Border EMPTY_BORDER = new EmptyBorder(1, 1, 1, 1);
	static private final Border BLACK_BORDER = LineBorder.createBlackLineBorder();
	static private final Border RED_BORDER = new LineBorder(Color.RED);

	private final Format format;
	private final JLabel renderComponent;
	private final JTextField editorComponent;
	private volatile transient Class<?> columnClass;


	/**
	 * Creates a new instance. Note that this implementation uses a subclass of JPanel as render component, which deactivates
	 * some standard method implementations for performance reasons.
	 * @param format the format converting the table cell values into text and vice versa, or {@code null}
	 */
	public TableCellFormater (final Format format) {
		this.format = format;
		this.editorComponent = new JTextField();
		this.renderComponent = new JLabel() {
			static private final long serialVersionUID = 1L;


			public void invalidate () {}


			public void validate () {}


			public void revalidate () {}


			public void repaint () {}
		};
	}


	/**
	 * Returns the format text for the given value.
	 * @param value the value or {@code null}
	 * @return the format text
	 */
	private String formatValue (final Object value) {
		if (this.format == null) {
			return Objects.toString(value, "");
		} else {
			try {
				return this.format.format(value);
			} catch (final Exception exception) {
				return "";
			}
		}
	}


	/**
	 * {@inheritDoc}
	 */
	public Component getTableCellRendererComponent (final JTable table, final Object value, final boolean selected, final boolean focused, final int rowIndex, final int columnIndex) {
		this.columnClass = this.format == null ? String.class : table.getColumnClass(columnIndex);

		final Border border = focused ? UIManager.getBorder(selected ? "Table.focusSelectedCellHighlightBorder" : "Table.focusCellHighlightBorder") : EMPTY_BORDER;
		final Color foreground = selected ? table.getSelectionForeground() : table.getForeground();
		final Color background = selected ? table.getSelectionBackground() : table.getBackground();

		return this.getCellRendererComponent(table, value, border, foreground, background, null);
	}


	/**
	 * {@inheritDoc}
	 */
	public Component getTreeCellRendererComponent (final JTree tree, final Object value, final boolean selected, final boolean expanded, final boolean leaf, final int row, final boolean focused) {
		this.columnClass = this.format == null ? String.class : (value == null ? Object.class : value.getClass());

		final Border border = new EmptyBorder(UIManager.getInsets("Tree.rendererMargins"));
		final Color foreground = selected ? UIManager.getColor("Tree.selectionForeground") : tree.getForeground();// use "Tree.textForeground"?
		final Color background = selected ? UIManager.getColor("Tree.selectionBackground") : tree.getBackground();// use "Tree.textBackground"?
		Icon icon = UIManager.getIcon(leaf ? "Tree.leafIcon" : (expanded ? "Tree.openIcon" : "Tree.closedIcon"));
		if (!tree.isEnabled()) UIManager.getLookAndFeel().getDisabledIcon(tree, icon);

		return this.getCellRendererComponent(tree, value, border, foreground, background, icon);
	}


	/**
	 * Primes the renderer component for use.
	 */
	private Component getCellRendererComponent (final JComponent parent, final Object value, final Border border, final Color foreground, final Color background, final Icon icon) {
		this.renderComponent.setText(this.formatValue(value));
		this.renderComponent.setOpaque(true);
		this.renderComponent.setBorder(border);
		this.renderComponent.setFont(parent.getFont());
		this.renderComponent.setForeground(foreground);
		this.renderComponent.setBackground(background);
		this.renderComponent.setComponentOrientation(parent.getComponentOrientation());
		this.renderComponent.setHorizontalAlignment(Number.class.isAssignableFrom(this.columnClass) ? SwingConstants.RIGHT : SwingConstants.LEFT);
		this.renderComponent.setEnabled(parent.isEnabled());
		this.renderComponent.setIcon(icon);
		return this.renderComponent;
	}


	/**
	 * {@inheritDoc}
	 */
	public Component getTableCellEditorComponent (final JTable table, final Object value, final boolean selected, final int rowIndex, final int columnIndex) {
		this.columnClass = this.format == null ? String.class : table.getColumnClass(columnIndex);
		return this.getCellEditorComponent(table, value);
	}


	/**
	 * {@inheritDoc}
	 */
	public Component getTreeCellEditorComponent (final JTree tree, final Object value, final boolean selected, final boolean expanded, final boolean leaf, final int rowIndex) {
		this.columnClass = this.format == null ? String.class : (value == null ? Object.class : value.getClass());
		return this.getCellEditorComponent(tree, value);
	}


	/**
	 * Primes the editor component for use.
	 */
	private Component getCellEditorComponent (final JComponent parent, final Object value) {
		this.editorComponent.setText(this.formatValue(value));
		this.editorComponent.setBorder(BLACK_BORDER);
		this.editorComponent.setFont(parent.getFont());
		this.editorComponent.setForeground(parent.getForeground());
		this.editorComponent.setBackground(parent.getBackground());
		this.editorComponent.setComponentOrientation(parent.getComponentOrientation());
		this.editorComponent.setHorizontalAlignment(Number.class.isAssignableFrom(this.columnClass) ? SwingConstants.RIGHT : SwingConstants.LEFT);
		return this.editorComponent;
	}


	/**
	 * {@inheritDoc}
	 * @throws IllegalStateException if cell validation fails
	 * @see #stopCellEditing()
	 */
	public Object getCellEditorValue () {
		final String text = this.editorComponent.getText();
		final Object value;
		if (this.format == null) {
			value = text;
		} else {
			try {
				value = this.format.parseObject(text);
			} catch (final Exception exception) {
				throw new IllegalStateException(exception);
			}
		}

		if (value == null | this.columnClass == null) return value;
		if (this.columnClass.isAssignableFrom(value.getClass())) return value;

		// make sure the correct number type is returned
		if (value instanceof Number & Number.class.isAssignableFrom(this.columnClass)) {
			try {
				@SuppressWarnings("unchecked")
				final Class<? extends Number> type = (Class<? extends Number>) this.columnClass;
				return Reflections.cast((Number) value, type);
			} catch (final IllegalArgumentException exception) {
				// do nothing
			}
		}
		return value;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean isCellEditable (final EventObject event) {
		if (!super.isCellEditable(event) | !(event instanceof MouseEvent)) return false;

		return ((MouseEvent) event).getClickCount() > 1;
	}


	/**
	 * {@inheritDoc}
	 */
	public boolean stopCellEditing () {
		if (this.format != null) {
			final String text = this.editorComponent.getText();
			try {
				this.format.parseObject(text);
			} catch (final Exception exception) {
				this.editorComponent.setBorder(RED_BORDER);
				return false;
			}
		}
		return super.stopCellEditing();
	}
}