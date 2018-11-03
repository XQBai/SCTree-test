# SCTree-test
SCTree test is an algorithm that can statistical detect the hidden structure of high-dimensional single-cell dataset, which the intrinsic structure may be linear structure of branched structure. Based on the tools of spiked matrix model and random matrix theory, SCTree construct the discordance matrix by transforming the distance between any pair of cells used Gromov-Farris transform. Then the discordance matrix can be formulated as the spiked matrix model. Therefore, we detect the spiked matrix model of discordance matrix whether signal or not. SCTree test have successfully identify the branch event of real RNA-seq data, meanwhile, by combining with the reconstruct method (such as Wishbone), it also can identify the multiple bifurcations of cellular process.  
### Installation and dependencies

1. SCTree test depends on a number of python3 packages available on pypi and these dependencies are listed in setup.py.  Most packages can be installed expect the special [Wishbone algorithm](https://github.com/ManuSetty/wishbone). Wishbone algorithm can be installed using:
```
   $> git clone git://github.com/ManuSetty/Palantir.git
   $> cd Palantir
   $> sudo -H pip3 install .   
```
    
2. SCTree test has been implemented in Python3. In code profile, SCTree.py can be downloaded to your own workspace. Then the function can be used by follow procedures:

```
    import SCTree
    data = sc.from_csv('~./data/linear.csv', data_type = 'sc-seq', normalize = False)
    p-value, SNR = sc.test(data, sampling_size, startcell, 'sc-seq')
    
```
If data need to be normalized, then the normlize = True, otherwise, False. "sampling_size" can be choiced as the number smaller than the data size. There need a given start cell "startcell". The result of "SNR" means the signal-to-noise ratio, the data has linear structure when "p-value" is greater than 0.05, otherwise, the data is branched structure.

