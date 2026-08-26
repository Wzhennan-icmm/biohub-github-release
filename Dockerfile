FROM rust:1.82-bookworm AS build
WORKDIR /src
COPY biohub-rs/Cargo.toml biohub-rs/Cargo.lock ./
COPY biohub-rs/src ./src
RUN cargo build --locked --release

FROM debian:bookworm-slim
LABEL org.opencontainers.image.title="BioHub"
LABEL org.opencontainers.image.description="Reproducible command-line utilities for plant genomics"
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates r-base samtools mafft pal2nal \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/biohub /usr/local/bin/biohub
COPY r /usr/local/share/biohub/r
COPY recipes /usr/local/share/biohub/recipes
COPY examples /usr/local/share/biohub/examples
COPY LICENSE NOTICE /usr/local/share/doc/biohub/
ENTRYPOINT ["biohub"]
CMD ["--help"]
