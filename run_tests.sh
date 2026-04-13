#!/bin/sh
source ./.venv/bin/activate
pytest tests/
cd go
make test
cd ..

