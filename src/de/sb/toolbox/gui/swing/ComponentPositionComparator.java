package de.sb.toolbox.gui.swing;

import java.awt.Component;
import java.awt.geom.Point2D;
import java.util.Comparator;
import java.util.Objects;
import de.sb.toolbox.Copyright;


/**
 * AWT/Swing component comparator based on screen location for use with a sorting focus traversal policy, allowing focus
 * traversal to ignore container hierarchy.
 * @see javax.swing.SortingFocusTraversalPolicy
 */
@Copyright(year = 2013, holders = "Sascha Baumeister")
public class ComponentPositionComparator implements Comparator<Component> {

	/**
	 * First compares the vertical positioning of both components; the ones placed higher up on the screen will be traversed
	 * before those displayed lower. In case both are only horizontally differing, the ones displayed more to the left of the
	 * screen will be traversed before those displayed more to the right. In case components are overlaying, their names are
	 * used to decide which ones are traversed first.
	 * @param leftComponent the left component, or {@code null}
	 * @param rightComponent the right component, or {@code null}
	 * @return {@code 1} if the left component shall be traversed before the right, {@code -1}if the right component shall be
	 *         traversed before the left, or {@code 0}if both components are considered equal.
	 */
	public int compare (final Component leftComponent, final Component rightComponent) {
		if (leftComponent == null & rightComponent == null) return 0;
		if (leftComponent == null) return -1;
		if (rightComponent == null) return 1;

		final Point2D leftOrigin = leftComponent.getLocationOnScreen();
		final Point2D rightOrigin = rightComponent.getLocationOnScreen();
		if (leftOrigin.getY() > rightOrigin.getY()) return 1;
		if (leftOrigin.getY() < rightOrigin.getY()) return -1;
		if (leftOrigin.getX() > rightOrigin.getX()) return 1;
		if (leftOrigin.getX() < rightOrigin.getX()) return -1;
		return Objects.toString(leftComponent.getName()).compareTo(Objects.toString(rightComponent.getName()));
	}
}