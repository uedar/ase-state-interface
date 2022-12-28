docker build \
    --build-arg pseudo_dir='./gncpp_pack' \
    --build-arg state_src='./state-5.6.10.tgz' \
    --no-cache \
    .