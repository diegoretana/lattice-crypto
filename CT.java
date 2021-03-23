package ABB;

import java.io.Serializable;

public class CT implements Serializable
{
    public long c0;
    public LVector c1;
    
    public CT(long c0, LVector c1)
    {
        this.c0 = c0;
        this.c1 = c1;
    }    
}
