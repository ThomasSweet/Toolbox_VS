package de.sb.toolbox.gui.swing;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.Arrays;
import java.util.Objects;
import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import de.sb.toolbox.Copyright;


/**
 * A vector area is a region that displays a vector value using a horizontal lollipop representation. It can be set to allow
 * editing or to be read-only. It's mode determines the editing semantics.
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class VectorArea extends JPanel {
	static private final long serialVersionUID = 1L;



	/**
	 * A vector area's editing tool.
	 */
	static public enum Tool {

		/**
		 * Clicking or dragging the mouse alters a single vector element.
		 */
		SLIDER,

		/**
		 * Clicking or dragging the mouse alters possibly multiple vector elements.
		 */
		PENCIL,

		/**
		 * Clicking or dragging the mouse negates possibly multiple vector elements.
		 */
		MIRROR,

		/**
		 * Clicking or dragging the mouse resets possibly multiple vector elements to zero.
		 */
		ERASER,

		/**
		 * Dragging the mouse to the left or right shifts all vector elements.
		 */
		SHIFTER
	}


	private volatile Tool editTool;
	private volatile Color scaleColor;
	private volatile float verticalScaleFactor;
	private volatile int highlightModule;
	private volatile float[] valueVector;


	/**
	 * Creates a new instance.
	 */
	public VectorArea () {
		this.editTool = Tool.SLIDER;
		this.verticalScaleFactor = 0.9f;
		this.valueVector = new float[1];

		final MouseEventListener listener = new MouseEventListener();
		this.addMouseListener(listener);
		this.addMouseMotionListener(listener);
		this.setPreferredSize(new Dimension(128, 64));
	}


	/**
	 * Fills all values with zeroes, and refreshes the draw panel.
	 */
	public void reset () {
		Arrays.fill(this.valueVector, 0);
		this.repaint();
	}


	/**
	 * Returns a copy of the value vector.
	 * @return the value vector
	 */
	public float[] getValueVector () {
		return this.valueVector.clone();
	}


	/**
	 * Sets the value vector by copying the elements of the given vector.
	 * @param valueVector the value vector
	 * @throws NullPointerException if the given argument is {@code null}
	 * @throws IllegalArgumentException if the given argument is empty
	 */
	public void setValueVector (final float[] valueVector) {
		if (valueVector.length == 0) throw new IllegalArgumentException();
		if (Arrays.equals(valueVector, this.valueVector)) return;

		if (this.valueVector.length != valueVector.length) {
			this.valueVector = new float[valueVector.length];
		}
		System.arraycopy(valueVector, 0, this.valueVector, 0, valueVector.length);
		this.repaint();

		this.fireStateChanged();
	}


	/**
	 * Returns the edit tool.
	 * @return the edit tool
	 */
	public Tool getEditTool () {
		return this.editTool;
	}


	/**
	 * Sets the edit tool.
	 * @param editTool the edit tool
	 * @throws NullPointerException if the given argument is {@code null}
	 */
	public void setEditTool (final Tool editTool) {
		if (editTool != this.editTool) {
			if (editTool == null) throw new NullPointerException();
			final Tool oldEditMode = this.editTool;
			this.editTool = editTool;
			this.firePropertyChange("editTool", oldEditMode, editTool);
		}
	}


	/**
	 * Returns the horizontal scale factor, which is the reciprocal of the value vector's length.
	 * @return the horizontal scale factor
	 */
	protected float getHorizontalScaleFactor () {
		return 1f / this.valueVector.length;
	}


	/**
	 * Returns the vertical scale factor, which is the reciprocal of the maximum visible element value.
	 * @return the vertical scale factor
	 */
	public float getVerticalScaleFactor () {
		return this.verticalScaleFactor;
	}


	/**
	 * Sets the vertical scale factor, which is the maximum visible element value.
	 * @param verticalScaleFactor the vertical scale factor
	 */
	public void setVerticalScaleFactor (final float verticalScaleFactor) {
		if (verticalScaleFactor != this.verticalScaleFactor) {
			final float oldVerticalScaleFactor = this.verticalScaleFactor;
			this.verticalScaleFactor = verticalScaleFactor;
			this.repaint();
			this.firePropertyChange("verticalScaleFactor", oldVerticalScaleFactor, verticalScaleFactor);
		}
	}


	/**
	 * Returns the color used to display the coordinate system. Unless this color is set explicitly, it is derived from this
	 * component's background color.
	 * @return the vertical scale color
	 */
	public Color getScaleColor () {
		if (this.scaleColor != null) return this.scaleColor;

		final float[] hsb = Color.RGBtoHSB(this.getBackground().getRed(), this.getBackground().getGreen(), this.getBackground().getBlue(), null);
		final int rgb = Color.HSBtoRGB(hsb[0], 0, hsb[2] > 0.5f ? 0.5f * (1 - hsb[2]) : 0.5f * (1 + hsb[2]));
		return new Color(rgb);
	}


	/**
	 * Sets the color used to display the coordinate system. Unless this color is set explicitly, it is derived from this
	 * component's background color.
	 * @return scaleColor the vertical scale color
	 */
	public void setScaleColor (final Color scaleColor) {
		if (!Objects.equals(scaleColor, this.scaleColor)) {
			final Color oldVerticalScaleColor = this.scaleColor;
			this.scaleColor = scaleColor;
			this.repaint();
			this.firePropertyChange("verticalScaleColor", oldVerticalScaleColor, scaleColor);
		}
	}


	/**
	 * Returns the highlight module. Every value element whose index fulfills the condition<br />
	 * {@code index % highlightModulo == 0} will be highlighted; a zero value means element highlighting is off.
	 * @return the highlight module
	 */
	public int getHighlightModule () {
		return this.highlightModule;
	}


	/**
	 * Sets the highlight module. Every value element whose index that fulfills the condition<br />
	 * {@code index % highlightModulo == 0} will be highlighted; a zero value means element highlighting is off.
	 * @param highlightModule the highlight module
	 */
	public void setHighlightModule (final int highlightModule) {
		if (highlightModule != this.highlightModule) {
			final int oldHighlightModule = this.highlightModule;
			this.highlightModule = highlightModule;
			this.repaint();
			this.firePropertyChange("highlightModule", oldHighlightModule, highlightModule);
		}
	}


	/**
	 * Adds a value change listener to this vector area.
	 * @param listener the change listener to be added
	 */
	public void addChangeListener (final ChangeListener listener) {
		this.listenerList.add(ChangeListener.class, listener);
	}


	/**
	 * Removes a value change listener from this vector area.
	 * @param listener the change listener to be removed
	 */
	public void removeChangeListener (final ChangeListener listener) {
		this.listenerList.remove(ChangeListener.class, listener);
	}


	/**
	 * Returns an array of all the value change listeners.
	 * @return the change listeners, or an empty array for none
	 */
	public ChangeListener[] getChangeListeners () {
		return this.listenerList.getListeners(ChangeListener.class);
	}


	/**
	 * Send a change event, whose source is this vector area, to all change listeners that have registered interest in change
	 * events.
	 */
	protected void fireStateChanged () {
		final ChangeEvent changeEvent = new ChangeEvent(this);
		for (final ChangeListener listener : this.getChangeListeners()) {
			listener.stateChanged(changeEvent);
		}
	}


	/**
	 * {@inheritDoc}
	 */
	protected void paintComponent (final Graphics pen) {
		super.paintComponent(pen);

		for (int index = 0; index < this.valueVector.length; ++index) {
			this.paintComponent(pen, index, false);
		}
	}


	/**
	 * Paints the horizontal component with the given index using the given graphics. If {@code clear} is {@code true}, then the
	 * component is first wiped by filling it with the background color. Then negative and positive horizontal unit axis are
	 * added in scale color. Finally, a horizontal zero axis and a lollipop for the vector element with the given index are
	 * added in foreground color.
	 * @param pen the graphics object to protect
	 * @param index the component index
	 * @param clear whether or not to clear the component before painting
	 */
	protected void paintComponent (final Graphics pen, final int index, final boolean clear) {
		final int halfHeight = this.getHeight() >> 1;
		final float floatUnitX = this.getWidth() * this.getHorizontalScaleFactor();
		final float floatUnitY = halfHeight * this.getVerticalScaleFactor();
		final int unitX = Math.round(floatUnitX);
		final int unitY = Math.round(floatUnitY);

		final Graphics clone = pen.create();
		try {
			clone.translate(Math.round(floatUnitX * index), halfHeight);
			if (clear) {
				clone.setColor(this.getBackground());
				clone.fillRect(0, -halfHeight, unitX, this.getHeight());
			}
			clone.setColor(this.getScaleColor());
			clone.drawLine(0, -unitY, unitX, -unitY);
			clone.drawLine(0, +unitY, unitX, +unitY);
			clone.setColor(this.getForeground());
			clone.drawLine(0, 0, unitX, 0);

			final int x = unitX >> 1;
			final int y = Math.round(-floatUnitY * this.valueVector[index]);
			final int r = Math.round(0.25f * floatUnitX);
			if (r < Math.abs(y)) {
				clone.drawLine(x, 0, x, y + (y < 0 ? +r : -r));
			}
			clone.drawOval(x - r, y - r, (r << 1), (r << 1));
			if (this.highlightModule == 0 || index % this.highlightModule != 0) {
				clone.fillOval(x - r, y - r, (r << 1), (r << 1));
			}
		} finally {
			clone.dispose();
		}
	}


	/**
	 * Modifies the value vector by setting the given element.
	 * @param elementIndex the element index
	 * @param elementValue the element value
	 * @throws IllegalArgumentException if the given element index is out of bounds
	 */
	protected void setValueElement (final int elementIndex, final float elementValue) {
		if (elementIndex < 0 | elementIndex >= this.valueVector.length) throw new IllegalArgumentException();
		if (elementValue == this.valueVector[elementIndex]) return;

		this.valueVector[elementIndex] = elementValue;
		this.paintComponent(this.getGraphics(), elementIndex, true);
		this.fireStateChanged();
	}


	/**
	 * Modifies the value vector by negating the given element.
	 * @param elementIndex the element index
	 * @throws IllegalArgumentException if the given element index is out of bounds
	 */
	protected void negateValueElement (final int elementIndex) {
		if (elementIndex < 0 | elementIndex >= this.valueVector.length) throw new IllegalArgumentException();
		if (this.valueVector[elementIndex] == 0) return;

		this.valueVector[elementIndex] = -this.valueVector[elementIndex];
		this.paintComponent(this.getGraphics(), elementIndex, true);
		this.fireStateChanged();
	}


	/**
	 * Modifies the value vector by zeroing the given element.
	 * @param elementIndex the element index
	 * @throws IllegalArgumentException if the given element index is out of bounds
	 */
	protected void zeroValueElement (final int elementIndex) {
		if (elementIndex < 0 | elementIndex >= this.valueVector.length) throw new IllegalArgumentException();
		if (this.valueVector[elementIndex] == 0) return;

		this.valueVector[elementIndex] = 0;
		this.paintComponent(this.getGraphics(), elementIndex, true);
		this.fireStateChanged();
	}


	/**
	 * Modifies the value vector by shifting all it's elements by the given amount to the right.
	 * @param shift the shift amount
	 */
	protected void shiftValueElements (int shift) {
		final int halfLength = this.valueVector.length >> 1;
		while (shift < -halfLength)
			shift += this.valueVector.length;
		while (shift >= halfLength)
			shift -= this.valueVector.length;
		if (shift == 0) return;

		final float[] swap = new float[Math.abs(shift)];
		if (shift < 0) {
			System.arraycopy(this.valueVector, 0, swap, 0, swap.length);
			System.arraycopy(this.valueVector, swap.length, this.valueVector, 0, this.valueVector.length - swap.length);
			System.arraycopy(swap, 0, this.valueVector, this.valueVector.length - swap.length, swap.length);
		} else {
			System.arraycopy(this.valueVector, this.valueVector.length - swap.length, swap, 0, swap.length);
			System.arraycopy(this.valueVector, 0, this.valueVector, swap.length, this.valueVector.length - swap.length);
			System.arraycopy(swap, 0, this.valueVector, 0, swap.length);
		}

		final Graphics pen = this.getGraphics();
		for (int index = 0; index < this.valueVector.length; ++index) {
			this.paintComponent(pen, index, true);
		}
		this.fireStateChanged();
	}



	/**
	 * Reacts on the associated parent's mouse events to provide editing and tool-tip functionality.
	 */
	private class MouseEventListener implements MouseListener, MouseMotionListener {
		static private final String TOOLTIP_FORMAT = "[%d -> %.3f]";

		private volatile int dragIndex = -1;


		/**
		 * {@inheritDoc}
		 */
		public void mousePressed (final MouseEvent event) {
			this.mousePressedOrDragged(event.getX(), event.getY());
		}


		/**
		 * {@inheritDoc}
		 */
		public void mouseDragged (final MouseEvent event) {
			this.mousePressedOrDragged(event.getX(), event.getY());
		}


		/**
		 * {@inheritDoc}
		 */
		public void mouseReleased (final MouseEvent event) {
			this.mouseReleasedOrExited();
		}


		/**
		 * {@inheritDoc}
		 */
		public void mouseExited (final MouseEvent event) {
			this.mouseReleasedOrExited();
		}


		/**
		 * {@inheritDoc}
		 */
		public void mouseMoved (final MouseEvent event) {
			final int segmentIndex = this.segmentIndex(event.getX());
			final String toolTipText = segmentIndex >= 0 & segmentIndex < VectorArea.this.valueVector.length ? String.format(TOOLTIP_FORMAT, segmentIndex, VectorArea.this.valueVector[segmentIndex]) : null;
			VectorArea.this.setToolTipText(toolTipText);
		}


		/**
		 * {@inheritDoc}
		 */
		public void mouseClicked (final MouseEvent event) {}


		/**
		 * {@inheritDoc}
		 */
		public void mouseEntered (final MouseEvent event) {}


		/**
		 * Alter the vector value and it's rendering. The operation returns immediately if the parent component is not enabled,
		 * of if any of the given coordinates is out of bounds.
		 * @param x the horizontal event coordinate
		 * @param y the vertical event coordinate
		 */
		private void mousePressedOrDragged (final int x, final int y) {
			final int width = VectorArea.this.getWidth();
			final int height = VectorArea.this.getHeight();
			if (x < 0 | x >= width | y < 0 | y >= height | !VectorArea.this.isEnabled()) return;

			final int segmentIndex = this.segmentIndex(x);
			switch (VectorArea.this.editTool) {
				case SLIDER: {
					final float elementValue = (1 - 2 * y / (float) height) / VectorArea.this.verticalScaleFactor;
					if (this.dragIndex == -1) this.dragIndex = segmentIndex;
					VectorArea.this.setValueElement(this.dragIndex, elementValue);
					break;
				}
				case PENCIL: {
					final float elementValue = (1 - 2 * y / (float) height) / VectorArea.this.verticalScaleFactor;
					VectorArea.this.setValueElement(segmentIndex, elementValue);
					this.dragIndex = segmentIndex;
					break;
				}
				case MIRROR: {
					if (segmentIndex != this.dragIndex) VectorArea.this.negateValueElement(segmentIndex);
					this.dragIndex = segmentIndex;
					break;
				}
				case ERASER: {
					VectorArea.this.zeroValueElement(segmentIndex);
					break;
				}
				default: {
					if (this.dragIndex != -1) VectorArea.this.shiftValueElements(segmentIndex - this.dragIndex);
					this.dragIndex = segmentIndex;
					break;
				}
			}
		}


		/**
		 * Reset mouse handling.
		 * @param event the mouse event
		 */
		private void mouseReleasedOrExited () {
			this.dragIndex = -1;
			VectorArea.this.setToolTipText(null);
		}


		/**
		 * Returns the index of the segment that is associated with the given x-coordinate.
		 * @param x the horizontal coordinate
		 * @return the associated segment index
		 */
		private int segmentIndex (final int x) {
			final int segmentCount = VectorArea.this.valueVector.length;
			final int index = x * segmentCount / VectorArea.this.getWidth();
			return Math.max(0, Math.min(segmentCount - 1, index));
		}
	}
}