package de.sb.toolbox.math;

public class PrimitivePolar {

	// TODO: use for addition!
	static public float[] add (final float abs1, final float arg1, final float abs2, final float arg2) {
		final float abs = (float) Math.sqrt(abs1 * abs1 + abs2 * abs2 + 2 * abs1 * abs2 * (float) Math.cos(arg1 - arg2));
		final float arg = arg1 - (float) Math.acos((abs1 * abs1 - abs2 * abs2 + abs * abs) / (2 * abs1 * abs));
		return new float[] { abs, arg };
	}
}
