package ABB;


import java.security.SecureRandom;
import com.wolfram.jlink.*;
import java.math.BigDecimal;
import java.math.RoundingMode;

public class Sample
{
    private static int n = SystemParameters.n;
    private static int q = SystemParameters.q;
    private static LVector t = null;
    
    private static SecureRandom sr = null;
    
    private static int scale = SystemParameters.scale;
    
    private static double t()
    {
        return Math.log(n)/Math.log(2);
    }
    
    public static int SampleZ(double c, double s)
    {
        long min = (long)Math.ceil( c - s * t());
        long max = (long)Math.floor(c + s * t());
        int x=0;

        KernelLink ml = SystemParameters.ml;
        
        String expression;
        try
        {
            expression = "min = " + String.valueOf(min) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "max = " + String.valueOf(max) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "c = " + String.valueOf(c) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "s = " + String.valueOf(s) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "NestWhile["
                    + "x = RandomInteger[{min,max}];"
                    + "rho = Exp[-Pi*(x-c)^2/s^2],"
                    + "x = RandomInteger[{min,max}];"
                    + "rho = Exp[-Pi*(x-c)^2/s^2],"
                    + "!( (rho > 0) && (rho <= 1) ) ];"
                    + "x";
            ml.evaluate(expression);
            ml.waitForAnswer();
            x = ml.getInteger();

        }catch(MathLinkException e)
        {
            System.out.println("MathLinkException occurred : " + e.getMessage());
        }
        
        return x;
    }
    
    public static LVector SampleD(LMatrix b, BigDecimal s, LVector c) 
    {
        int i, zi;
        BigDecimal cp, sp;
        LVector ci;
        LVector vi;
        
        ci = c;
        int dim = b.basis[0].getDim();
        vi = new LVector(dim); // zero vector
        
        for(i= dim-1; i>=0; i--)
        {
            cp = ( b.gso[i].punto(ci) ).divide( b.gso[i].Norm2, scale, RoundingMode.HALF_UP );
            sp = s.divide( new BigDecimal( b.gso[i].Norm ), scale, RoundingMode.HALF_UP);
            zi = SampleZ(cp.doubleValue(),sp.doubleValue() );
            
            ci = ci.sub( b.basis[i].mult(zi) );
            vi = vi.add( b.basis[i].mult(zi) );
        }

        return vi;
    }
    
    public static LVector Sample_Zm1(int m1, BigDecimal sigma)
    {
        long[][]  Im1 = new long[m1][m1];
        long[] zerom1 = new long[m1];
        
        LMatrix I = new LMatrix(Im1);
        I.setIdentity();
        LVector zero = new LVector(zerom1);
        
        return SampleD(I,sigma,zero);
    }
    
    public static LVector SamplePre(LMatrix A, LMatrix Ta, LVector u, BigDecimal sigma)
    {
        int q = SystemParameters.q;
        KernelLink ml = SystemParameters.ml;
        long[] tarray;
            
        String expression;

        try
        {
            //ml.discardAnswer();  // enable when execution is not automatic
            expression = "A = " + A.getString() + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "u = " + u.getString() + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "t = Table[Subscript[x, i], {i, " + String.valueOf(Ta.rows) + "}];";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "rules = FindInstance[A.t == u, t, Modulus -> " + String.valueOf(q) + "];"
                       + "Flatten[ t /. rules ]";
            ml.evaluate(expression);
            ml.waitForAnswer();

            tarray = ml.getLongArray1();

            t = new LVector(tarray);

            System.out.println(" t : ");
            t.print(); System.out.println("");

        }catch(MathLinkException e)
        {
            System.out.println("MathLinkException occurred : " + e.getMessage());
        }
        
        LVector v = SampleD(Ta, sigma, t.mult(-1));
        
        return t.add(v);
    }
    
    public static LVector SampleLeft2(PublicParam PP, LMatrix M1, LMatrix Ta, LVector target, LVector e2)
    {
        LVector e1 = SamplePre2( Ta, PP.sigma, target);
        
        return e1.augment(e2);
    }
    
    private static LVector SamplePre2( LMatrix Ta, BigDecimal sigma, LVector target )
    {
        LVector v = SampleD( Ta, sigma, target.mult(-1) );
        
        return target.add( v );
    }
    

    public static LVector SampleLeft(LMatrix A, LMatrix M1, LMatrix Ta, LVector u, BigDecimal sigma)
    {
        //int m  = SystemParameters.m;
        int m1 = M1.columns;
        
        
        
        System.out.println("sample_zm1");
        LVector e2 = Sample_Zm1(m1,sigma);
        //System.out.println(" e2 : ");
        //e2.print(); System.out.println("");
        //System.out.println(" m1 : " + m1 + "   sigma : " + sigma);
        //System.out.println("|| e2 || : " + e2.norm() );System.out.println("\n");
        
        //long[] iy = {5037669, 4816524, 4329788};
        LVector y  = u.sub(M1.Mult(e2));
        y.reduceModule(q);
        //System.out.println(" y : ");
        //y.print(); System.out.println("\n");
        
        System.out.println(" samplePre ");
        LVector e1 = SamplePre(A, Ta, y, sigma);
        //System.out.println(" e1 : ");
        //e1.print(); System.out.println("\n");
        //System.out.println(" || e1 || : " + e1.norm());
        
        return e1.augment(e2);
    }
    
}