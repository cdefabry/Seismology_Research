/*	Grav2d.java

	Usage:	java grav2d tol infile outfile

	Parameter:	tol=0.1		mGal tolerance to match data to

	Compile:	javac grav2d.java

	Input File Format:

		nel niter drho(g/cc) meter?
		x grav(mGal)
		.
		.
		. (nel data lines)
	
	If meter?==true x is in meters; if anything else or not present
	x is in km.
*/

import java.io.*;
import java.util.StringTokenizer;

class grav2d
	{
	static DataInputStream dis;
	static StringTokenizer st;
	@SuppressWarnings("deprecation")
	public static void main(String[] argv)
		{
		double tol=0.1;
		int nel, niter, i, it=1, l, m, ma;
		double drho, grav=0.0, g=6.67e-8, t1=1e8;
		double dp, dg, ba, x, x1, x2, th;
		double max, x0;
		boolean meter;
		String s;
		if (argv.length != 3)
			{
			System.out.println("java grav2d tol infile outfile");
			return;
			}
		try
			{
			tol = new Double(argv[0]).doubleValue();
			}
		catch (NumberFormatException nfe)
			{
			System.out.println("grav2d: improper tol=" + argv[0] + ", abort.");
			return;
			}
		FileInputStream fis;
		try
			{
			fis = new FileInputStream(argv[1]);
			}
		catch (FileNotFoundException fnfe)
			{
			System.out.println("grav2d: can't open input file "
				+ argv[1] + ", abort.");
			return;
			}
		dis = new DataInputStream(fis);
		
		/* scanf("%d %d %lf", &nel, &niter, &drho); */
		try
			{
			st = new StringTokenizer(dis.readLine());
			nel = Integer.parseInt(st.nextToken());
			niter = Integer.parseInt(st.nextToken());
			drho = new Double(st.nextToken()).doubleValue();
			}
		catch (IOException ioe)
			{
			System.out.println("grav2d: couldn't read first line from input file "
				+ argv[1] + ", abort.");
			return;
			}
		catch (NumberFormatException nfe)
			{
			System.out.println("grav2d: couldn't parse nel, niter, drho from first line of input file "
				+ argv[1] + ", abort.");
			return;
			}
		if (st.hasMoreTokens())
			{
			try
				{
				meter = new Boolean(st.nextToken()).booleanValue();
				}
			catch (NumberFormatException nfe)
				{
				meter = false;
				}
			}
		else
			meter = false;
		Strip el[] = new Strip[nel];
		for (i=0; i<nel; i++)
			{
			el[i] = new Strip();
			try
				{
				st = new StringTokenizer(dis.readLine());
				el[i].x = new Double(st.nextToken()).doubleValue();
				el[i].gmeas = new Double(st.nextToken()).doubleValue();
				}
			catch (IOException ioe)
				{
				System.out.println("grav2d: couldn't read line "
					+ (i+2) + " from input file " + argv[1] + ", abort.");
				return;
				}
			catch (NumberFormatException nfe)
				{
				System.out.println("grav2d: couldn't parse line "
					+ (i+2) + " for x, gmeas from input file "
					+ argv[1] + ", abort.");
				return;
				}
			if (meter)
				el[i].x /= 1000;
			}
		/* prepare x1, x2, ganom as in gravin.c */
		x0 = el[0].x;
		/*el[0].x = 0.0;*/
		max = el[0].gmeas;
		for (i=1; i<nel; i++)
			{
			/*el[i].x -= x0;*/
			/* for drho<0, normalize ganom to maximum gmeas */
			if (drho < 0)
				{
				if (el[i].gmeas > max)
					max = el[i].gmeas;
				}
			/* for drho>0, normalize ganom to minimum gmeas */
			else
				{
				if (el[i].gmeas < max)
					max = el[i].gmeas;
				}
			}
		el[0].x2 = (el[0].x + el[1].x)/2.0 - x0;
		el[0].x1 = -el[0].x2;
		el[nel-1].x1 = (el[nel-2].x + el[nel-1].x)/2.0 - x0;
		el[nel-1].x2 = el[nel-1].x + (el[nel-1].x - el[nel-1].x1 - x0) - x0;
		for (i=1; i<nel-1; i++) 
			{
			el[i].x2 = (el[i].x + el[i+1].x)/2.0 - x0;
			el[i].x1 = (el[i-1].x + el[i].x)/2.0 - x0;
			}
		for (i=0; i<nel; i++)
			{
			el[i].b = el[i].x2 - el[i].x1;
			if (drho < 0)
				el[i].ganom = el[i].gmeas - max - 1.0;
			else
				el[i].ganom = el[i].gmeas - max + 1.0;
			}
		FileOutputStream fos;
		try
			{
			fos = new FileOutputStream(argv[2]);
			}
		catch (IOException ioe)
			{
			System.out.println("grav2d: can't create output file "
				+ argv[2] + ", abort.");
			return;
			}
		PrintStream ps = new PrintStream(fos);
		ps.println("" + nel + " Elements, " + niter + " Iterations Max, "
			+ tol + " Tolerance, drho=" + drho);
		ps.println("");
				
		do
			{
			if (it > 1)
				{
				for (m=0; m<nel; m++)
					{
					dp = el[m].dep;
					ba = el[m].b;
					x = ba/2.0;
					grav += tdim(ba, x, drho, dp);
					for (ma=0; ma<m; ma++)
						{
						x1 = el[m].x1 + el[m].b/2.0 - el[ma].x1;
						x1 = x1>0 ? x1 : -x1;
						x2 = el[m].x1 + el[m].b/2.0 - el[ma].x2;
						x2 = x2>0 ? x2 : -x2;
						x = x2>x1 ? x2 : x1;
						dp = el[ma].dep;
						ba = el[ma].b;
						grav += tdim(ba, x, drho, dp);
						}
					el[m].gcalc = grav;
					grav = 0.0;
					}
				}
			else
				{
				for (i=0; i<nel; i++)
					{
					el[i].dep = 0.0;
					el[i].gcalc = 0.0;
					}
				}
	
			ps.println("Iteration " + it + ":");
			ps.println("Distance\tGmeas\tDepth\tX1\t\tX2\t\tGanom\tGcalc");
			for (i=0; i<nel; i++)
				ps.println(getSignif("" + el[i].x, 5) + "\t"
				+ getSignif("" + el[i].gmeas) + "\t"
				+ getSignif("" + el[i].dep) + "\t"
				+ getSignif("" + el[i].x1) + "\t"
				+ getSignif("" + el[i].x2) + "\t"
				+ getSignif("" + el[i].ganom) + "\t"
				+ getSignif("" + el[i].gcalc));
			ps.println("");
			System.out.println("grav2d: finished iteration " + it);
	
			if (it == niter)
				{
				ps.println("Convergence not reached by max. iteration "
					+ it);
				System.out.println("Convergence not reached by max. iteration "
					+ it);
				return;
				}
	
			l = 0;
			for (i=0; i<nel; i++)
				{
				dg = el[i].ganom - el[i].gcalc;
				if ((dg>0?dg:-dg) <= tol)
					{
					l++;
					if (l == nel)
						{
						ps.println("Convergence at iteration " + it);
						System.out.println("Convergence at iteration " + it);
						return;
						}
					}
				th = dg/6.28318/g/drho/t1;
				el[i].dep += th;
				if (el[i].dep <= 0.0)
					el[i].dep = 0.0;
				}
			it++;
			}
		while (it <= niter);
		}
	
	static double tdim(double b, double x, double drho, double dp)
		{
		double t1, t2, g, x1, x2, x3, x4, anom;
		t1 = 1e8;
		t2 = 1e-4;
		if (dp <= t2)
			anom = 0.0;
		else
			{
			g = 6.67e-8;
			x1 = x * Math.log(Math.sqrt((dp*dp + x*x)/(x*x)));
			x2 = (x-b) * Math.log(Math.sqrt((dp*dp+(x-b)*(x-b))/((x-b)*(x-b))));
			x3 = Math.atan2(x, dp);
			x4 = Math.atan2(x-b, dp);
			anom = 2.0*g*drho*(x1 - x2 + dp*(x3-x4))*t1;
			}
		return anom;
		}
	
	/* Method helping display of float numbers */
	static String getSignif(String flt)
		{
		return getSignif(flt, 4);
		}
	static String getSignif(String flt, int sig)
		{
		int i, dot, exp;
		float val, aval;
		if (sig < 8)
			sig = 8;
		if (flt.length() < sig)
			return flt;
		dot = flt.indexOf('.');
		exp = flt.indexOf('E');
		if (exp < 0 || exp >= flt.length())
			exp = flt.indexOf('e');
		if (dot < 1 || dot >= flt.length())
			return flt;
		try {val = Float.valueOf(flt).floatValue(); }
		catch (NumberFormatException nfe)
			{return flt; }
		aval = val < 0F ? -val : val;
		if (aval < 1F && (exp < 0 || exp > flt.length()))
			{
			if (val < 0F)
				return flt.substring(0, sig);
			return flt.substring(0, sig-1);
			}
		if (val < 0F)
			{
			if (dot > (sig-1))
				return flt.substring(0, dot-1);
			if (dot == (sig-1))
				return flt.substring(0, sig);
			if (dot == 2 && exp > dot && exp < flt.length())
				return flt.substring(0, sig-2) + flt.substring(exp);
			return flt.substring(0, sig-1);
			}
		if (dot > (sig-2))
			return flt.substring(0, dot-1);
		if (dot == (sig-2))
			return flt.substring(0, sig-1);
		if (dot == 1 && exp > dot && exp < flt.length())
			return flt.substring(0, sig-3) + flt.substring(exp);
		return flt.substring(0, sig-2);
		}
	} /* end of grav2d class */

class Strip
	{
	double x;
	double x1;
	double x2;
	double gmeas;
	double ganom;
	double dep;
	double b;
	double gcalc;
	}
