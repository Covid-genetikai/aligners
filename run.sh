docker run \
    -v $(pwd)/tree/src:/src \
    -v $(pwd)/datasets/ncbi:/src/data \
    -it bio:latest sh