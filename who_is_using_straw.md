# Software using Straw for `.hic` files

## Direct dependency / direct call

  ---------------------------------------------------------------------------------------------------------------------------
  Software          Purpose               Straw usage                    Link
  ----------------- --------------------- ------------------------------ ----------------------------------------------------
  Mustache          Chromatin loop        `hic-straw`; directly reads    https://github.com/ay-lab/mustache
                    calling               `.hic`                         

  Peakachu          ML chromatin loop     `hic-straw` dependency for     https://github.com/tariks/peakachu
                    calling               `.hic` input                   

  SELFISH           Differential          `hic-straw` dependency         https://github.com/ay-lab/selfish
                    chromatin                                            
                    interactions                                         

  diffDomain        Differential TAD      `hic-straw` dependency         https://github.com/Tian-Dechao/diffDomain
                    detection                                            

  HiCLift           Hi-C genome-assembly  `hic-straw` for `.hic` input   https://github.com/XiaoTaoWang/HiCLift
                    liftover                                             

  HIPPS-DIMES       3D genome             `hic-straw` for `.hic` input   https://github.com/anyuzx/HIPPS-DIMES
                    reconstruction                                       

  HiCRep            Hi-C reproducibility  `strawr` used by `hic2mat()`   https://github.com/TaoYang-dev/hicrep

  CALDER / CALDER2  Domains, hierarchy,   Imports `strawr`; accepts      https://github.com/CSOgroup/CALDER2
                    subcompartments       `.hic`                         

  STRIDE /          Hi-C similarity       Requires `hic-straw` for       https://pypi.org/project/hicstride/
  hicstride         metric                `.hic`                         

  scHiCTools        Single-cell Hi-C      Adapts Straw for `.hic`        https://github.com/liu-bioinfo-lab/scHiCTools
                    analysis              reading                        

  PEKORA            Hi-C analysis         `hic-straw` for `.hic` input   https://github.com/sXperfect/pekora

  EPInformer        Enhancer-promoter     Environment includes           https://github.com/pinellolab/EPInformer
                    prediction            `hicstraw`                     

  HiCExperiment     Bioconductor Hi-C     Imports `strawr`               https://bioconductor.org/packages/HiCExperiment
                    infrastructure                                       

  HicAggR           Aggregate Hi-C        Imports `strawr`               https://bioconductor.org/packages/HicAggR
                    analysis                                             

  mariner           Paired genomic        Imports `strawr`               https://bioconductor.org/packages/mariner
                    interaction analysis                                 

  plotgardener      Genomic visualization Imports `strawr`               https://bioconductor.org/packages/plotgardener

  trackViewer       Genomic track         Imports `strawr`               https://bioconductor.org/packages/trackViewer
                    visualization                                        

  AutoHiC           Automated chromosome  Uses `hic-straw` to extract    https://doi.org/10.1093/nar/gkae749
                    assembly              `.hic` matrices                

  4DGBWorkflow /    Dynamic 3D genome     Uses `hic-straw` to parse      https://pmc.ncbi.nlm.nih.gov/articles/PMC12338655/
  4DHiC             modeling              `.hic`                         

  pC-SAC            3D genome             Calls                          https://pmc.ncbi.nlm.nih.gov/articles/PMC11995266/
                    reconstruction        `hicstraw.getMatrixZoomData`   

  ProDy chromatin   Structural/dynamics   Includes Straw-based `.hic`    https://github.com/prody/ProDy
                    analysis              reader                         
  ---------------------------------------------------------------------------------------------------------------------------

## Straw as documented preprocessing

  ---------------------------------------------------------------------------------------------------------
  Software /        Purpose           Straw             Link
  workflow                            preprocessing     
  ----------------- ----------------- ----------------- ---------------------------------------------------
  TADCompare        Differential TAD  Official workflow https://github.com/dozmorovlab/TADCompare
                    analysis          uses Straw to     
                                      extract `.hic`    
                                      matrices before   
                                      analysis          

  SpectralTAD       TAD calling       Official vignette https://bioconductor.org/packages/SpectralTAD
                                      uses Straw to     
                                      convert `.hic` to 
                                      matrix input      

  PHi-C2            Dynamic 3D genome Paper extracts    https://doi.org/10.1093/bioinformatics/btac626
                    modeling          contact matrices  
                                      from `.hic` with  
                                      Straw             

  TopDom workflows  TAD calling       Published         https://pmc.ncbi.nlm.nih.gov/articles/PMC9768916/
                                      workflows use     
                                      Straw `.hic`      
                                      extraction before 
                                      TopDom            

  GENOVA workflows  Hi-C analysis and Published         https://github.com/robinweide/GENOVA
                    visualization     analyses dump     
                                      `.hic` contacts   
                                      with Straw before 
                                      GENOVA            

  TimeCompare /     Temporal TAD      Published         https://doi.org/10.1038/s41597-022-01508-x
  TADCompare        analysis          workflow extracts 
  workflows                           `.hic` matrices   
                                      with Straw before 
                                      comparison        
  ---------------------------------------------------------------------------------------------------------

## Notes

-   **Direct dependency / direct call** means the software itself
    imports, requires, incorporates, or explicitly calls Straw,
    `hic-straw`/`hicstraw`, or `strawr`.
-   **Documented preprocessing** means Straw is used upstream to extract
    matrices from `.hic`, while the downstream software consumes the
    resulting matrix rather than necessarily calling Straw itself.
