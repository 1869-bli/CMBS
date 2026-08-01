"""CMBS GUI: type content, preview the colour code, export PNG, decode images back.

Requires Pillow.  Run with:  python gui.py
"""

import os
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox

from PIL import Image, ImageTk

from cmbs import (
    __version__,
    build_grid,
    decode_image,
    encode_to_image,
    payload_capacity,
    COLOR_NAMES,
    COLOR_GLYPHS,
    PALETTE,
    GRID_SIZE,
    ECC_LEVELS,
    DecodeError,
)

try:
    import cv2
    from cmbs.camera import decode_frame

    _HAS_CAMERA = True
except ImportError:
    cv2 = None
    _HAS_CAMERA = False

PREVIEW_SIZE = 420


class CMBSApp:
    def __init__(self, root):
        self.root = root
        root.title("CMBS v%s - Colour Manifested Byte Storage" % __version__)
        try:
            icon = Image.open(os.path.join(os.path.dirname(__file__), "cmbslogo.png"))
            root.iconphoto(True, ImageTk.PhotoImage(icon.resize((64, 64))))
        except Exception:
            pass

        self.level = tk.StringVar(value="M")
        self.cell = tk.IntVar(value=20)
        self.border = tk.IntVar(value=4)
        self.last_grid = None

        self._build_widgets()
        self._update_capacity()

    # -- layout ------------------------------------------------------------

    def _build_widgets(self):
        outer = tk.Frame(self.root, padx=10, pady=10)
        outer.pack(fill="both", expand=True)

        left = tk.Frame(outer)
        left.pack(side="left", fill="both", expand=True)

        tk.Label(left, text="Content (UTF-8):", anchor="w").pack(fill="x")
        self.text = tk.Text(left, width=42, height=14, wrap="word", font=("Consolas", 10))
        self.text.pack(fill="both", expand=True)
        self.text.bind("<KeyRelease>", lambda e: self._update_capacity())

        controls = tk.Frame(left)
        controls.pack(fill="x", pady=(8, 0))

        tk.Label(controls, text="Error correction:").pack(side="left")
        menu = tk.OptionMenu(controls, self.level, *ECC_LEVELS, command=lambda v: self._update_capacity())
        menu.pack(side="left", padx=(4, 12))

        tk.Label(controls, text="Cell px:").pack(side="left")
        tk.Spinbox(controls, from_=8, to=60, textvariable=self.cell, width=4).pack(side="left", padx=(2, 12))
        tk.Label(controls, text="Border:").pack(side="left")
        tk.Spinbox(controls, from_=1, to=10, textvariable=self.border, width=4).pack(side="left", padx=(2, 12))

        self.cap_label = tk.Label(controls, text="")
        self.cap_label.pack(side="right")

        buttons = tk.Frame(left)
        buttons.pack(fill="x", pady=(8, 0))
        tk.Button(buttons, text="Generate preview", command=self.generate).pack(side="left")
        tk.Button(buttons, text="Export PNG", command=self.export).pack(side="left", padx=6)
        tk.Button(buttons, text="Decode image...", command=self.decode).pack(side="left")
        tk.Button(buttons, text="Scan with camera", command=self.scan_camera).pack(side="left", padx=(6, 0))
        self.status = tk.Label(left, text="", anchor="w", fg="#444")
        self.status.pack(fill="x", pady=(6, 0))

        right = tk.Frame(outer, width=PREVIEW_SIZE + 16, height=PREVIEW_SIZE + 16)
        right.pack(side="left", padx=(12, 0), fill="y")
        right.pack_propagate(False)

        self.canvas = tk.Canvas(right, width=PREVIEW_SIZE, height=PREVIEW_SIZE, bg="white",
                                highlightthickness=1, highlightbackground="#bbb")
        self.canvas.pack()
        self.canvas.create_text(PREVIEW_SIZE // 2, PREVIEW_SIZE // 2, text="Generate a preview",
                                fill="#999", font=("Segoe UI", 12))

        legend = tk.Frame(outer)
        legend.pack(side="left", padx=(12, 0), fill="y")
        tk.Label(legend, text="Palette:", font=("Segoe UI", 9, "bold")).pack(anchor="w")
        for i, name in enumerate(COLOR_NAMES):
            row = tk.Frame(legend)
            row.pack(anchor="w", pady=1)
            tk.Label(row, text=COLOR_GLYPHS[i], font=("Segoe UI Emoji", 12)).pack(side="left")
            tk.Label(row, text="  %03d = %s" % (i, name), font=("Segoe UI", 9), anchor="w").pack(side="left")

    # -- behaviour ----------------------------------------------------------

    def _payload(self):
        return self.text.get("1.0", "end-1c").encode("utf-8")

    def _update_capacity(self, *_):
        cap = payload_capacity(self.level.get())
        used = len(self._payload())
        ok = used <= cap
        color = "#060" if ok else "#c00"
        self.cap_label.config(text="%d / %d bytes" % (used, cap), fg=color)

    def generate(self, *_):
        data = self._payload()
        try:
            grid = build_grid(data, self.level.get())
        except ValueError as exc:
            messagebox.showerror("CMBS", str(exc))
            return
        self.last_grid = grid
        img = encode_to_image(data, self.level.get(), cell=self.cell.get(), border=self.border.get())
        thumb = img.resize((PREVIEW_SIZE, PREVIEW_SIZE), Image.NEAREST)
        self._photo = ImageTk.PhotoImage(thumb)
        self.canvas.delete("all")
        self.canvas.create_image(PREVIEW_SIZE // 2, PREVIEW_SIZE // 2, image=self._photo)
        self.status.config(text="%d bytes encoded (level %s)" % (len(data), self.level.get()))

    def export(self):
        if self.last_grid is None:
            messagebox.showinfo("CMBS", "Generate a preview first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png", filetypes=[("PNG image", "*.png")],
            initialfile="cmbs.png")
        if not path:
            return
        encode_to_image(self._payload(), self.level.get(),
                        cell=self.cell.get(), border=self.border.get()).save(path)
        self.status.config(text="exported %s" % path)

    def scan_camera(self):
        if not _HAS_CAMERA:
            messagebox.showinfo(
                "CMBS",
                "Camera scanning needs OpenCV. Install it with:\n\n"
                "    pip install opencv-python")
            return
        CameraWindow(self)

    def decode(self):
        path = filedialog.askopenfilename(
            filetypes=[("PNG image", "*.png"), ("Image", "*.png *.jpg *.jpeg *.bmp")])
        if not path:
            return
        try:
            data, rot = decode_image(path)
        except DecodeError as exc:
            messagebox.showerror("CMBS", str(exc))
            self.status.config(text="decode failed")
            return
        try:
            text = data.decode("utf-8")
            if not text.isprintable() or "\x00" in text:
                raise UnicodeDecodeError("utf-8", data, 0, len(data), "not printable")
        except UnicodeDecodeError:
            text = data.hex()
        self.text.delete("1.0", "end")
        self.text.insert("1.0", text)
        self.status.config(text="decoded %d bytes (rotation %d)" % (len(data), rot))
        self._update_capacity()


class CameraWindow:
    """Live webcam preview that continuously tries to decode a CMBS code."""

    PREVIEW_W, PREVIEW_H = 640, 480

    def __init__(self, app):
        self.app = app
        self.running = True
        self.queue = queue.Queue()

        self.top = tk.Toplevel(app.root)
        self.top.title("Scan CMBS with camera")
        self.top.protocol("WM_DELETE_WINDOW", self.stop)
        self.canvas = tk.Canvas(self.top, width=self.PREVIEW_W, height=self.PREVIEW_H,
                                bg="black", highlightthickness=0)
        self.canvas.pack(padx=8, pady=(8, 0))
        self.status = tk.Label(
            self.top, text="Hold the code flat, well lit, and roughly level.",
            fg="#444", font=("Segoe UI", 10))
        self.status.pack(pady=4)
        tk.Button(self.top, text="Stop", command=self.stop).pack(pady=(0, 8))

        self._photo = None
        self.thread = threading.Thread(target=self._capture, daemon=True)
        self.thread.start()
        self.top.after(40, self._poll)

    def stop(self):
        self.running = False
        self.top.destroy()

    def _capture(self):
        cap = cv2.VideoCapture(0, cv2.CAP_DSHOW)
        if not cap.isOpened():
            cap = cv2.VideoCapture(0)
        if not cap.isOpened():
            self.queue.put(("error", "No camera found."))
            return
        frame_no = 0
        while self.running:
            ok, frame = cap.read()
            if not ok:
                continue
            frame = cv2.resize(frame, (self.PREVIEW_W, self.PREVIEW_H),
                               interpolation=cv2.INTER_AREA)
            result = None
            if frame_no % 4 == 0:
                try:
                    result = decode_frame(frame)
                except Exception:
                    result = None
            frame_no += 1
            self.queue.put(("frame", frame, result))
        cap.release()

    def _poll(self):
        if not self.running:
            return
        try:
            item = self.queue.get_nowait()
        except queue.Empty:
            pass
        else:
            if item[0] == "error":
                messagebox.showerror("CMBS", item[1])
                self.stop()
                return
            _, frame, result = item
            rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
            self._photo = ImageTk.PhotoImage(Image.fromarray(rgb))
            self.canvas.delete("all")
            self.canvas.create_image(self.PREVIEW_W // 2, self.PREVIEW_H // 2, image=self._photo)
            if result is not None:
                data, rot = result
                try:
                    text = data.decode("utf-8")
                    if not text.isprintable() or "\x00" in text:
                        raise UnicodeDecodeError("utf-8", data, 0, len(data), "not printable")
                except UnicodeDecodeError:
                    text = data.hex()
                self.app.text.delete("1.0", "end")
                self.app.text.insert("1.0", text)
                self.app._update_capacity()
                self.app.status.config(text="scanned %d bytes from camera (rotation %d)" % (len(data), rot))
                self.status.config(text="Decoded %d bytes!" % len(data))
                self.top.after(400, self.stop)
                return
            self.status.config(text="Scanning...")
        self.top.after(40, self._poll)


def main():
    root = tk.Tk()
    try:
        root.tk.call("tk", "scaling", 1.0)
    except tk.TclError:
        pass
    CMBSApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
