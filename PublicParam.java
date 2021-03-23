package ABB;

import java.math.BigDecimal;
import java.util.concurrent.ThreadLocalRandom;
import java.io.Serializable;

public class PublicParam implements Serializable
{
    public LMatrix A0;
    public LMatrix A1;
    public LMatrix B;
    public LVector[] u;
    public int n;
    public int q;
    public int m;
    public String irreduciblePolynomial;
    public double alpha;
    public BigDecimal sigma;
    public LVector[] target;
    public LVector[] e2;
    public LVector id;
    
    private int i, j; // counter
    private int me = 24; // multiencrypt block
    
    
    
    public PublicParam(LMatrix a0, LMatrix a1, LMatrix b, 
            int n, int m, int q, String poly, double alpha, BigDecimal sigma)
    {
        A0 = a0;
        A1 = a1;
        B = b;
        this.n = n;
        this.q = q;
        this.m = m;
        this.irreduciblePolynomial = poly;
        this.alpha = alpha;
        this.sigma = sigma;
        
        this.id = createID();
        
        this.target  = computeSet_target();
        this.e2 = computeSet_e2();
        this.u  = computeSet_u();
    }
    
    private LVector createID()
    {
        long[] iid = new long[ n ];
        
        for(i=0; i<n; i++)
            iid[i] = 1;
        
        return new LVector(iid);
    }
    
    private LVector[] computeSet_target()
    {
        LVector[] sety = new LVector[me];
        long[] it = new long[m];
        
        for(j=0; j<me; j++)
        {
            for(i=0; i<m; i++)
                it[i] = ThreadLocalRandom.current().nextInt(0,2);
            sety[j] = new LVector(it);
            //System.out.print(" t " + j + " : ");
            //sety[j].print();
            //System.out.println("");
        }

        return sety;
    }
    
    private LVector[] computeSet_e2()
    {
        LVector[] sete2 = new LVector[me];
        
        for(j=0; j<me; j++)
        {
            sete2[j] = Sample.Sample_Zm1(m, sigma.divide( new BigDecimal("10") ) );
            //System.out.println("e2 " + j + " : ");
            //sete2[j].print();
        }
        
        return sete2;
    }
    
    private LVector[] computeSet_u()
    {
        SystemParameters.irreduciblePoly = this.irreduciblePolynomial;
        SystemParameters.n = this.n;
        SystemParameters.q = this.q;
        
        LMatrix HidB = Scheme.H(id, 1);
        HidB = HidB.ModMult(B, q);
        HidB = A1.addition(HidB);
        
        LVector[] setu = new LVector[me];
        LVector y;
        
        for(j=0; j<me; j++)
        {
            y = A0.ModMult( target[j] );
            setu[j] = y.add( HidB.Mult(e2[j]) );
            setu[j].reduceModule(q);
            //System.out.println(" u " + j + " : ");
            //setu[j].print();
        }
        
        return setu;
    }
    
}
