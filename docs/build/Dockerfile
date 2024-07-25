FROM python:3.11-alpine AS builder

COPY ./requirements.txt .
RUN pip install -r requirements.txt
COPY ./sphinx/ ./src/
COPY ./module/ ./src/module/
RUN sphinx-build -b html ./src ./src/_build/html

FROM nginx:1.21.3
## WARNING: Without a reverse proxy, the subpath of /usr/share/nginx/html
## must match the ingress path in the deployment. 
COPY --from=builder /src/_build/html/ /usr/share/nginx/html/docs
