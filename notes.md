# Installation on udd machine

Installing pytorch

```
pip install \
    --extra-index-url=https://pypi.nvidia.com \
    "cudf-cu13==25.10.*" "dask-cudf-cu13==25.10.*" "cuml-cu13==25.10.*" \
    "cugraph-cu13==25.10.*" "nx-cugraph-cu13==25.10.*" "cuxfilter-cu13==25.10.*" \
    "cucim-cu13==25.10.*" "pylibraft-cu13==25.10.*" "raft-dask-cu13==25.10.*" \
    "cuvs-cu13==25.10.*" "nx-cugraph-cu13==25.10.*"
```


```
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu130
```

Installing transformers

```
pip install transformers
```

```
conda install nvidia::libcurand-dev
conda install nvidia/label/cuda-13.0.0::cuda-toolkit
```
