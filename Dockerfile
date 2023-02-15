FROM alpine:latest

RUN apk add --no-cache \
  bash \
  perl \
  make \
  g++ \
  bzip2-dev \
  zlib-dev

WORKDIR /app
COPY . . 
RUN ./install_krakenuniq.sh /usr/local/bin

CMD ["bash"]
