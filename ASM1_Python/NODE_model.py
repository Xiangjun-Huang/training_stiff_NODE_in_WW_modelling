import torch.nn as nn

### Define NODE model ###

class neuralODE(nn.Module):

    def __init__(self, nf, m, s, dm, ds, mx, mn, dmx, dmn, yts):
        super(neuralODE, self).__init__()
        self.nf = nf
        self.ymu = m
        self.ystd = s
        self.dymu = dm
        self.dystd = ds
        self.ymax = mx
        self.ymin = mn
        self.dymax = dmx
        self.dymin = dmn
        self.ytscale = yts

        self.net = nn.Sequential(
            nn.Linear(self.nf, 50),
            # nn.Tanh(),
            nn.GELU(),
            nn.Linear(50, 50),
            # nn.Tanh(),
            nn.GELU(),
            nn.Linear(50, 50),
            # nn.Tanh(),
            nn.GELU(),
            nn.Linear(50, self.nf))

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, mean=0, std=0.5)

    def forward(self, t, y):

        # choose the methods: direct, standardisation, max-min normalisation, or equation scaling
        # direct output, without normalisation
        # y = self.net(y)

        # z-score standardisation
        y = (y - self.ymu) / self.ystd   
        y = self.net(y)
        y = y * self.dystd + self.dymu

        # max-min normalisation
        # ymax_ymin = self.ymax - self.ymin
        # y[:, ymax_ymin == 0.0] = 0.0
        # y[:, ymax_ymin != 0.0] = (y[:, ymax_ymin != 0.0] - self.ymin[ymax_ymin != 0]) / ymax_ymin[ymax_ymin != 0]
        # y = self.net(y)
        # y = y * (self.dymax - self.dymin) + self.dymin

        # equation scaling, refer to Kim:"stiff neural ordinary equations"
        # loss function must be matched with dividend : /yscale
        # y = self.net(y) 
        # y = y * self.ytscale

        return y