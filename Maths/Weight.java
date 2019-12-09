package NumCalc;

import java.util.*;

		/**  Weighting strategy 
		*/
		public enum Weight
		{
			/**  No weights (constant) 
			*/
			Constant,

			/**  Linear increasing weights 
			*/
			Linear,

			/**  Polynomial increasing weights (quadratic, cubic...) 
			*/
			Polynomial;

			public int getValue()
			{
				return this.ordinal();
			}

			public static Weight forValue(int value)
			{
				return values()[value];
			}
		}
		