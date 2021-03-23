package ABB;

import com.wolfram.jlink.*;
import java.security.SecureRandom;
import java.util.Date;
import java.security.NoSuchAlgorithmException;
import java.util.concurrent.ThreadLocalRandom;


public class RandomGen 
{   
    private static int q = SystemParameters.q;
    private static int m = SystemParameters.m;
    private static int n = SystemParameters.n;
    
    public static long Gaussian_Zq( )
    {
        KernelLink ml = SystemParameters.ml;
        long x=0;
        
        String expression, stralpha;
        try
        {
            stralpha = Double.toString(SystemParameters.alpha);
            stralpha = stralpha.replace("E", "*10^");
            
            expression = "alpha = " + stralpha + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "sd = alpha/Sqrt[2*Pi];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "X = Mod[ RandomVariate[ NormalDistribution[0,sd] ], 1 ];";
            ml.evaluate(expression);
            ml.discardAnswer();
            
            expression = "q = " + String.valueOf(SystemParameters.q) + ";";
            ml.evaluate(expression);
            ml.discardAnswer();

            expression = "Mod[Round[ q*X ],q]";
            ml.evaluate(expression);
            ml.waitForAnswer();
            x = ml.getLongInteger();

        }catch(MathLinkException e)
        {
            System.err.println(" GaussianZq MathLinkException occurred : " + e.getMessage());
        }
        
        return x;
    }
    
    public static LVector Gaussian_Zq_m( )
    {
        long[] v = new long[m];
        
        int i;
        for(i=0; i<m; i++)
            v[i] = Gaussian_Zq( );
        
        return new LVector(v);
    }
    
    public static LVector Random_Zq_n()
    {
        long[] v = new long[n];
        
        SecureRandom sr = null;
        
        try
        {
            sr = SecureRandom.getInstance("SHA1PRNG");
            sr.setSeed(new Date().getTime());
        }catch(NoSuchAlgorithmException nsae){}
        
        int i;
        for(i=0; i<n; i++)
            v[i] = sr.nextInt(q);
        
        return new LVector(v);
    }
    
     public static LMatrix UniformRandom_Zq_nxm(int n, int m)
    {
        long[][] array = new long[n][m];
        
        for(int i=0; i<n; i++)
            for(int j=0; j<m; j++)
                array[i][j] = ThreadLocalRandom.current().nextInt(0,q-1);

        return new LMatrix(array);
    }
    
    public static LMatrix Random_R()
    {
        long[][] R  = new long[m][m];
        
        SecureRandom sran = new SecureRandom();
        sran.setSeed(new Date().getTime());
        
        int i, j, randomNum;
        
        for(i=0; i<m; i++)
            for(j=0; j<m; j++)
            {
                randomNum = sran.nextInt(2);
                if(randomNum == 0)
                    R[i][j] = -1;
                else
                    R[i][j] = 1;
            }

        return new LMatrix(R);
    }   
}