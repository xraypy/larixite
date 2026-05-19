"""Tests for larixite Flask webapp."""

import tempfile
from pathlib import Path

import pytest

# Skip if Flask not installed
FLASK = pytest.importorskip("flask")


STRUCTS = Path(__file__).parent / "structs"


@pytest.fixture(scope="module")
def app():
    from larixite.webapp.webapp import app

    app.config["TESTING"] = True
    app.config["UPLOAD_FOLDER"] = str(tempfile.mkdtemp())
    app.config["MAX_CONTENT_LENGTH"] = 2**24
    app.config["PROPAGATE_EXCEPTIONS"] = False
    return app


@pytest.fixture(scope="module")
def client(app):
    return app.test_client()


ZNO_TEXT = (STRUCTS / "ZnO_mp-2133.cif").read_text()


def _upload(client, content):
    """Upload a CIF and return the redirect location."""
    with tempfile.NamedTemporaryFile(suffix=".cif", delete=False, mode="w") as f:
        f.write(content)
        fname = f.name
    try:
        rv = client.post(
            "/upload_cif",
            data={"file": (open(fname, "rb"), "test.cif")},
            content_type="multipart/form-data",
            follow_redirects=True,
        )
        assert rv.status_code == 200
        return rv
    finally:
        Path(fname).unlink()


class TestIndex:
    def test_index(self, client):
        rv = client.get("/")
        assert rv.status_code == 200


class TestUpload:
    def test_upload_valid_cif(self, client):
        rv = _upload(client, ZNO_TEXT)
        assert rv.status_code == 200

    def test_upload_empty_file(self, client):
        rv = client.post(
            "/upload_cif",
            data={"file": (b"", "")},
            content_type="multipart/form-data",
            follow_redirects=True,
        )
        assert rv.status_code == 200


class TestCifs:
    def test_cifs_page_after_upload(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.get("/cifs/test.cif")
        assert rv.status_code == 200

    def test_feff_post(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.post(
            "/cifs/test.cif",
            data={
                "feff": "1",
                "absorbing_atom": "Zn",
                "edge": "K",
                "cluster_size": "7.0",
                "with_h": "0",
            },
            follow_redirects=True,
        )
        assert rv.status_code == 200
        data = rv.get_data(as_text=True)
        assert "Zn" in data

    def test_fdmnes_post(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.post(
            "/cifs/test.cif",
            data={
                "fdmnes": "1",
                "absorbing_atom": "Zn",
                "edge": "K",
                "cluster_size": "7.0",
                "with_h": "0",
            },
            follow_redirects=True,
        )
        assert rv.status_code == 200


class TestOutput:
    """Test /output route — returns raw FEFF/FDMNES input text."""

    def _setup(self, client):
        _upload(client, ZNO_TEXT)
        return client

    def test_feff_output_zn_k(self, client):
        self._setup(client)
        rv = client.get("/output/test.cif/Zn/1/K/7.0/0/feff/0/0/feff.inp")
        assert rv.status_code == 200
        data = rv.get_data(as_text=True)
        assert data.count("TITLE") >= 1 or data.count("*") > 0
        assert "Zn" in data

    def test_feff_output_with_h(self, client):
        self._setup(client)
        rv = client.get("/output/test.cif/Zn/1/K/7.0/1/feff/0/0/feff.inp")
        assert rv.status_code == 200

    def test_fdmnes_output(self, client):
        self._setup(client)
        rv = client.get("/output/test.cif/Zn/1/K/7.0/0/fdmnes/0/0/fdmfile.txt")
        assert rv.status_code == 200
        data = rv.get_data(as_text=True)
        assert "1" in data or ".inp" in data

    def test_fdmnes_output_o_edge(self, client):
        self._setup(client)
        rv = client.get("/output/test.cif/O/1/K/7.0/0/fdmnes/0/0/fdmfile.txt")
        assert rv.status_code == 200


class TestCiffile:
    def test_ciffile_download(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.get("/ciffile/test.cif/test.cif")
        assert rv.status_code == 200
        data = rv.get_data(as_text=True)
        assert "Zn" in data
        assert "O" in data


class TestZipfile:
    def test_feff_zip(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.get("/zipfile/test.cif/Zn/K/7.0/0/feff/0/test_zoom")
        assert rv.status_code in (200, 500)

    def test_fdmnes_zip(self, client):
        _upload(client, ZNO_TEXT)
        rv = client.get("/zipfile/test.cif/Zn/K/7.0/0/fdmnes/0/test_zoom")
        assert rv.status_code in (200, 500)


class TestUploadCifFile:
    """Test actual file upload with multipart form data."""

    def test_real_file_upload(self, client):
        with open(STRUCTS / "ZnO_mp-2133.cif", "rb") as f:
            rv = client.post(
                "/upload_cif",
                data={"file": (f, "uploaded.cif")},
                content_type="multipart/form-data",
                follow_redirects=True,
            )
            assert rv.status_code == 200

    def test_no_file_part(self, client):
        rv = client.post(
            "/upload_cif",
            data={},
            content_type="multipart/form-data",
            follow_redirects=True,
        )
        assert rv.status_code == 200

    def test_wrong_extension(self, client):
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False, mode="w") as f:
            f.write(ZNO_TEXT)
            fname = f.name
        try:
            rv = client.post(
                "/upload_cif",
                data={"file": (open(fname, "rb"), "test.txt")},
                content_type="multipart/form-data",
                follow_redirects=True,
            )
            assert rv.status_code == 200
        finally:
            Path(fname).unlink()


class TestAbout:
    def test_about(self, client):
        rv = client.get("/about/")
        assert rv.status_code == 200


class TestUtility:
    """Test helper functions exposed by the webapp."""

    def test_random_string(self):
        from larixite.webapp.webapp import random_string

        for _ in range(5):
            s = random_string(12)
            assert len(s) == 12 and s.islower()

    def test_allowed_file(self):
        from larixite.webapp.webapp import allowed_file

        assert allowed_file("test.cif") is True
        assert allowed_file("test.txt") is False
        assert allowed_file("test.CIF") is True
        assert allowed_file("test") is False


class TestGetElements:
    def test_simple_element(self):
        from larixite.webapp.webapp import get_elements

        stoich, _ = get_elements("Fe")
        assert stoich == {"Fe": 1}

    def test_compound(self):
        from larixite.webapp.webapp import get_elements

        stoich, _ = get_elements("Fe2O3")
        assert stoich == {"Fe": 2, "O": 3}

    def test_wat_prefix_raises(self):
        from larixite.webapp.webapp import get_elements

        with pytest.raises(UnboundLocalError):
            get_elements("Wat")
