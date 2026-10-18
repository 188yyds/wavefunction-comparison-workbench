import os, sys, json, copy, csv, base64, math, shutil, traceback, re, html
from datetime import datetime
from pathlib import Path
from decimal import Decimal
from model import *
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data'/'matplotlib'))
os.environ.setdefault('QT_API','pyqt5')
from PyQt5 import QtCore, QtGui, QtWidgets as Q
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import font_manager
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():
    name=font_manager.FontProperties(fname=str(font)).get_name()
    if hasattr(font_manager.fontManager,'addfont'): font_manager.fontManager.addfont(str(font))
    else: font_manager.fontManager.ttflist.append(font_manager.FontEntry(fname=str(font),name=name))
    matplotlib.rcParams['font.family']=name
matplotlib.rcParams['axes.unicode_minus']=False

STATUS={'computed':'完成','cached':'已复用缓存','failed':'失败','incompatible':'不适用',
        'point':'开始该点','point_done':'该点完成','calculating':'计算中','finished':'阶段结束'}
STYLE={'实线':'-','虚线':'--','点线':':','点划线':'-.'}

def button(text,fn,parent_layout):
    b=Q.QPushButton(text);b.clicked.connect(fn);parent_layout.addWidget(b);return b
def line(value): return Q.QLineEdit(str(value))
def safe_name(text): return re.sub(r'[<>:"/\\|?*\x00-\x1f]','_',text)[:100]

class ResultCanvas(FigureCanvasQTAgg):
    failed=QtCore.pyqtSignal(str)
    def draw(self,*args,**kwargs):
        try:super().draw(*args,**kwargs)
        except Exception:self.failed.emit(traceback.format_exc())

class DefinitionsDialog(Q.QDialog):
    def __init__(self,parent=None):
        super().__init__(parent);self.setWindowTitle('物理量计算说明');self.resize(870,650)
        layout=Q.QVBoxLayout(self);bar=Q.QHBoxLayout();layout.addLayout(bar)
        self.search=line('');self.search.setPlaceholderText('搜索物理量、公式或区间');bar.addWidget(self.search,1)
        button('查找下一处',self.find_next,bar);self.search.returnPressed.connect(self.find_next)
        self.browser=Q.QTextBrowser();self.browser.setOpenExternalLinks(False);layout.addWidget(self.browser,1)
        QtGui.QFontDatabase.addApplicationFont(font_manager.findfont('DejaVu Sans'))
        sections=[]
        for paragraph in (ROOT/'物理量定义.md').read_text(encoding='utf-8').split('\n\n'):
            escaped=html.escape(paragraph).replace('\n','<br>')
            for character,letter in [('ᵢ','i'),('ⱼ','j'),('ᵣ','r'),('ₐ','a'),('ₖ','k')]:escaped=escaped.replace(character,'<sub>'+letter+'</sub>')
            for character,letter in [('ᵇ','b'),('ⁱ','i'),('ᴸ','L')]:escaped=escaped.replace(character,'<sup>'+letter+'</sup>')
            if paragraph.startswith('# '):sections.append('<h1>'+escaped[2:]+'</h1>')
            elif paragraph.startswith('## '):sections.append('<h2>'+escaped[3:]+'</h2>')
            else:sections.append('<p>'+escaped+'</p>')
        self.browser.setHtml('<html><body style="font-size:12pt;line-height:150%">'+''.join(sections)+'</body></html>')
        button('关闭',self.close,layout)
    def find_next(self):
        term=self.search.text().strip()
        if term and not self.browser.find(term):
            self.browser.moveCursor(QtGui.QTextCursor.Start);self.browser.find(term)

class BatchImportDialog(Q.QDialog):
    def __init__(self,project,group,parent=None,start=None):
        super().__init__(parent);self.project=project;self.group=group
        self.setWindowTitle('批量导入波函数 · '+group['name']);self.resize(920,650)
        layout=Q.QVBoxLayout(self);bar=Q.QHBoxLayout();layout.addLayout(bar)
        bar.addWidget(Q.QLabel('从 U 开始填入空位'));self.start=Q.QComboBox();self.start.addItems(group_grid(project,group));bar.addWidget(self.start)
        if start in group_grid(project,group):self.start.setCurrentText(start)
        button('选择多个文件',self.choose,bar);button('添加文件夹内的波函数',self.choose_folder,bar)
        hint=Q.QLabel('按列表顺序对应 U；可拖动文件行或用上移 / 下移调整。已有绑定跳过，超出空位时整批拒绝。')
        hint.setWordWrap(True);layout.addWidget(hint)
        self.files=Q.QListWidget();self.files.setSelectionMode(Q.QAbstractItemView.ExtendedSelection)
        self.files.setDragDropMode(Q.QAbstractItemView.InternalMove);layout.addWidget(self.files,1)
        controls=Q.QHBoxLayout();layout.addLayout(controls)
        button('上移',lambda:self.move_file(-1),controls);button('下移',lambda:self.move_file(1),controls)
        button('移除所选',self.remove_selected,controls);controls.addStretch()
        self.preview=Q.QTableWidget(0,2);self.preview.setHorizontalHeaderLabels(['对应 U','完整文件路径'])
        self.preview.setEditTriggers(Q.QAbstractItemView.NoEditTriggers);self.preview.horizontalHeader().setStretchLastSection(True)
        self.preview.setMaximumHeight(180);layout.addWidget(self.preview)
        self.summary=Q.QLabel();layout.addWidget(self.summary)
        buttons=Q.QDialogButtonBox(Q.QDialogButtonBox.Ok|Q.QDialogButtonBox.Cancel);layout.addWidget(buttons)
        self.import_button=buttons.button(Q.QDialogButtonBox.Ok);self.import_button.setText('按此顺序导入')
        buttons.button(Q.QDialogButtonBox.Cancel).setText('取消')
        buttons.accepted.connect(self.accept);buttons.rejected.connect(self.reject)
        self.files.model().rowsMoved.connect(self.refresh_preview);self.start.currentIndexChanged.connect(self.refresh_preview)
        self.refresh_preview()
    def paths(self):return [self.files.item(i).text() for i in range(self.files.count())]
    def add_paths(self,paths):
        for path in paths:
            item=Q.QListWidgetItem(os.path.abspath(path));item.setToolTip(item.text());self.files.addItem(item)
        self.refresh_preview()
    def choose(self):
        paths,_=Q.QFileDialog.getOpenFileNames(self,'批量选择波函数；下一步核对顺序','','波函数 (*.h5 *.hdf5)')
        self.add_paths(paths)
    def choose_folder(self):
        directory=Q.QFileDialog.getExistingDirectory(self,'添加当前文件夹中的 .h5 / .hdf5（不递归）')
        if directory:
            try:
                with os.scandir(directory) as entries:
                    self.add_paths([e.path for e in entries if e.is_file() and Path(e.name).suffix.lower() in ('.h5','.hdf5')])
            except OSError as e:Q.QMessageBox.warning(self,'无法读取文件夹',str(e))
    def move_file(self,direction):
        rows=sorted(self.files.row(item) for item in self.files.selectedItems())
        if not rows or (direction<0 and rows[0]==0) or (direction>0 and rows[-1]==self.files.count()-1):return
        for row in (rows if direction<0 else reversed(rows)):
            item=self.files.takeItem(row);self.files.insertItem(row+direction,item);item.setSelected(True)
        self.refresh_preview()
    def remove_selected(self):
        for row in sorted((self.files.row(item) for item in self.files.selectedItems()),reverse=True):self.files.takeItem(row)
        self.refresh_preview()
    def refresh_preview(self,*args):
        us=group_grid(self.project,self.group);us=us[us.index(self.start.currentText()):]
        slots=[u for u in us if not self.group['files'].get(u)];paths=self.paths();self.preview.setRowCount(len(paths))
        for i,path in enumerate(paths):
            for j,value in enumerate([slots[i] if i<len(slots) else '超出空位',path]):
                item=Q.QTableWidgetItem(value);item.setToolTip(value);self.preview.setItem(i,j,item)
        valid=bool(paths) and len(paths)<=len(slots)
        self.import_button.setEnabled(valid);self.summary.setText('已选 %d 个文件；可用 %d 个空位%s'%(len(paths),len(slots),'；请减少文件或扩展 U 网格' if len(paths)>len(slots) else ''))

class FileTable(Q.QTableWidget):
    dropped=QtCore.pyqtSignal(list,int)
    swapped=QtCore.pyqtSignal(int,int)
    def __init__(self):
        super().__init__(0,4)
        self.setHorizontalHeaderLabels(['U','文件名','完整路径','状态'])
        self.setAcceptDrops(True);self.viewport().setAcceptDrops(True)
        self.setSelectionBehavior(Q.QAbstractItemView.SelectRows)
        self.setSelectionMode(Q.QAbstractItemView.ExtendedSelection)
        self.setEditTriggers(Q.QAbstractItemView.NoEditTriggers)
        self.setDragEnabled(True)
        self.horizontalHeader().setStretchLastSection(True)
        self.setColumnWidth(0,75);self.setColumnWidth(1,235);self.setColumnWidth(2,260)
    def dragEnterEvent(self,e):
        if e.mimeData().hasUrls() or e.mimeData().hasFormat('application/x-wc-row'): e.acceptProposedAction()
    def dragMoveEvent(self,e): self.dragEnterEvent(e)
    def dropEvent(self,e):
        row=self.rowAt(e.pos().y())
        if e.mimeData().hasFormat('application/x-wc-row'):
            if row>=0: self.swapped.emit(int(bytes(e.mimeData().data('application/x-wc-row')).decode()),row)
        else:
            paths=[u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()]
            if paths: self.dropped.emit(paths,row)
        e.acceptProposedAction()
    def startDrag(self,actions):
        row=self.currentRow()
        if row<0: return
        mime=QtCore.QMimeData();mime.setData('application/x-wc-row',str(row).encode('ascii'))
        drag=QtGui.QDrag(self);drag.setMimeData(mime);drag.exec_(QtCore.Qt.MoveAction)

class ParamsDialog(Q.QDialog):
    def __init__(self,params,parent=None,overrides=False,grid_values=None):
        super().__init__(parent);self.setWindowTitle('组内设置' if overrides else '物理量与计算设置')
        self.resize(570,700);layout=Q.QVBoxLayout(self);form=Q.QFormLayout()
        scroll=Q.QScrollArea();body=Q.QWidget();body.setLayout(form);scroll.setWidget(body);scroll.setWidgetResizable(True);layout.addWidget(scroll)
        self.override=Q.QCheckBox('为本组覆盖公共物理设置');self.override.setChecked(bool(overrides))
        if overrides is not False: form.addRow(self.override)
        self.inputs={}
        fields=[('t1','t1'),('ratio','ratio'),('dimer','dimer（1 或 2）'),('nbar','非局域参考密度'),
          ('L','预期 L（0 = 从 MPS 检查）'),('N','预期 N（-1 = 从 MPS 检查）'),
          ('dataset','MPS 数据集'),('bulk_left','体区间左端（0 = 自动）'),('bulk_right','体区间右端（0 = 自动）'),
          ('string_left','String 左端（0 = 自动）'),('string_right','String 右端（0 = 自动）'),
          ('cut','半链与纠缠谱切口（0 = 中心）'),('entropy_coarse','块熵粗投影阈值'),
          ('entropy_fine','块熵细投影阈值'),('memory_mb','块熵临时内存限额 / MB'),('blas_threads','BLAS 线程数')]
        for key,label in fields:
            w=line(params[key]);self.inputs[key]=w;form.addRow(label,w)
        self.custom=Q.QCheckBox('本组使用独立 U 网格');self.custom.setChecked(grid_values is not None)
        self.ginputs=[line(v) for v in (grid_values or ['1','2','0.5'])]
        if overrides is not False:
            form.addRow(self.custom)
            for label,w in zip(['U 起点','U 终点','U 步长'],self.ginputs):form.addRow(label,w)
            def enable_physics(checked):
                for w in self.inputs.values():w.setEnabled(checked)
            def enable_grid(checked):
                for w in self.ginputs:w.setEnabled(checked)
            self.override.toggled.connect(enable_physics);self.custom.toggled.connect(enable_grid)
            enable_physics(self.override.isChecked());enable_grid(self.custom.isChecked())
        help=Q.QLabel('区间两端同时为 0 时自动确定。块熵不自动放宽精度；高阶累积量自动进行数值稳定性检查。')
        help.setWordWrap(True);layout.addWidget(help)
        bb=Q.QDialogButtonBox(Q.QDialogButtonBox.Ok|Q.QDialogButtonBox.Cancel)
        bb.accepted.connect(self.check_accept);bb.rejected.connect(self.reject);layout.addWidget(bb)
    def values(self):
        out={}
        integer={'dimer','L','N','bulk_left','bulk_right','string_left','string_right','cut','memory_mb','blas_threads'}
        for key,w in self.inputs.items():
            out[key]=w.text().strip() if key=='dataset' else int(w.text()) if key in integer else float(w.text())
        if not out['dataset']: raise ValueError('数据集名称不能为空')
        if out['dimer'] not in (1,2): raise ValueError('dimer 必须是 1 或 2')
        if not all(math.isfinite(v) for v in out.values() if isinstance(v,(int,float))):raise ValueError('参数必须有限')
        if out['t1']<=0 or out['ratio']<0 or out['L']<0 or out['N']<-1: raise ValueError('t1 > 0，ratio ≥ 0，L ≥ 0，N ≥ -1')
        if out['memory_mb']<64 or out['blas_threads']<1: raise ValueError('内存限额至少 64 MB，线程数至少 1')
        if not 0<out['entropy_fine']<out['entropy_coarse']<1: raise ValueError('要求 0 < 细阈值 < 粗阈值 < 1')
        for a,b in [('bulk_left','bulk_right'),('string_left','string_right')]:
            if not ((out[a]==out[b]==0) or (0<out[a]<out[b])):raise ValueError('区间端点须同时为 0，或 0 < 左端 < 右端')
        if out['cut']<0: raise ValueError('切口不能为负')
        return out
    def check_accept(self):
        try:
            self.values()
            if self.custom.isChecked():grid(*(w.text() for w in self.ginputs))
            self.accept()
        except Exception as e:Q.QMessageBox.warning(self,'设置错误',str(e))

class MainWindow(Q.QMainWindow):
    def __init__(self,project=None):
        super().__init__()
        if font.exists():
            font_id=QtGui.QFontDatabase.addApplicationFont(str(font))
            families=QtGui.QFontDatabase.applicationFontFamilies(font_id)
            if families:Q.QApplication.instance().setFont(QtGui.QFont(families[0],9))
        self.setWindowTitle('波函数对比工作台');self.resize(1460,920)
        self.project=project or new_project();self.running=False;self.rows=[];self.run_dir=None
        self.snapshot=None;self.buffer='';self.logfile=None;self.checked_status={};self.phase='normal'
        self.runtime=dict(python=sys.executable,julia=shutil.which('julia') or 'julia',julia_project='')
        cfg=ROOT/'data'/'runtime.json'
        if cfg.exists():
            try:self.runtime.update(json.loads(cfg.read_text(encoding='utf-8')))
            except Exception:pass
        self.page_scroll=Q.QScrollArea();self.page_scroll.setWidgetResizable(True);self.page_scroll.setFrameShape(Q.QFrame.NoFrame)
        self.setCentralWidget(self.page_scroll);container=Q.QWidget();container.setMinimumSize(900,600)
        self.page_scroll.setWidget(container);outer=Q.QVBoxLayout(container)
        top=Q.QHBoxLayout();outer.addLayout(top)
        title=Q.QLabel('波函数对比工作台');title.setStyleSheet('font-size:22px;font-weight:600;color:#163451');top.addWidget(title);top.addStretch()
        for text,fn in [('新方案',self.new),('打开方案',self.open_project),('保存方案',self.save_project),('载入示例',self.load_example),('运行环境',self.runtime_dialog)]:button(text,fn,top)
        button('物理量计算说明',self.show_definitions,top)
        self.tabs=Q.QTabWidget();outer.addWidget(self.tabs)
        setup=Q.QWidget();self.setup=setup;self.tabs.addTab(setup,'设置与导入');sl=Q.QVBoxLayout(setup)
        bar=Q.QHBoxLayout();sl.addLayout(bar);self.grid_inputs=[line(v) for v in self.project['grid']]
        for label,w in zip(['U 起点','终点','步长'],self.grid_inputs):bar.addWidget(Q.QLabel(label));w.setMaximumWidth(105);bar.addWidget(w)
        button('生成 / 更新 U 网格',self.apply_global_grid,bar);button('公共物理设置',self.common_params,bar)
        bar.addStretch();bar.addWidget(Q.QLabel('参考组'));self.reference=Q.QComboBox();bar.addWidget(self.reference)
        self.reference.currentIndexChanged.connect(self.change_reference)
        splitter=Q.QSplitter();sl.addWidget(splitter)
        files=Q.QWidget();fl=Q.QVBoxLayout(files);splitter.addWidget(files)
        groupsbar=Q.QHBoxLayout();fl.addLayout(groupsbar);self.groups=Q.QComboBox();groupsbar.addWidget(self.groups,1)
        self.groups.currentIndexChanged.connect(self.refresh_table)
        button('添加组',self.add_group,groupsbar);button('组设置',self.group_dialog,groupsbar);button('删除组',self.remove_group,groupsbar)
        self.group_note=Q.QLabel();self.group_note.setTextFormat(QtCore.Qt.PlainText);self.group_note.setWordWrap(True);fl.addWidget(self.group_note)
        self.table=FileTable();self.table.dropped.connect(self.import_files);self.table.swapped.connect(self.swap_rows);fl.addWidget(self.table,1)
        fb=Q.QHBoxLayout();fl.addLayout(fb)
        for text,fn in [('选择文件',self.choose_files),('批量导入波函数',self.batch_import),('替换选中行',self.replace_selected),('交换两行',self.swap_selected),('移除 / 留空',self.clear_selected),('检查对应关系',lambda:self.check_mapping())]:button(text,fn,fb)
        self.pending=Q.QListWidget();self.pending.setSelectionMode(Q.QAbstractItemView.ExtendedSelection);self.pending.setMaximumHeight(85);fl.addWidget(Q.QLabel('待分配文件（选中后可填入空位）'));fl.addWidget(self.pending)
        pb=Q.QHBoxLayout();fl.addLayout(pb);button('分配待选文件',self.assign_pending,pb);button('移除待选记录',self.remove_pending,pb)
        selection=Q.QWidget();ol=Q.QVBoxLayout(selection);splitter.addWidget(selection)
        ob=Q.QHBoxLayout();ol.addLayout(ob);ob.addWidget(Q.QLabel('最高可选累积量阶数'))
        self.order=Q.QSpinBox();self.order.setRange(6,32);self.order.setValue(self.project['common']['max_order']);ob.addWidget(self.order)
        self.order.valueChanged.connect(self.change_order)
        self.tree=Q.QTreeWidget();self.tree.setHeaderLabels(['选择物理量']);ol.addWidget(self.tree,1)
        self.tree.itemChanged.connect(self.selection_changed)
        presets=Q.QHBoxLayout();ol.addLayout(presets)
        button('快速预设',lambda:self.select_preset('quick'),presets);button('Scaling 预设',lambda:self.select_preset('scaling'),presets)
        button('全选',lambda:self.select_preset('all'),presets);button('清空',lambda:self.select_preset('none'),presets)
        hint=Q.QLabel('仅计算勾选项及必要依赖。中心块熵较耗内存；超过六阶的累积量需要高精度检查。')
        hint.setWordWrap(True);hint.setStyleSheet('color:#617085');ol.addWidget(hint);splitter.setSizes([950,390])
        results=Q.QWidget();self.tabs.addTab(results,'结果与导出');rl=Q.QVBoxLayout(results)
        rb=Q.QHBoxLayout();rl.addLayout(rb);self.obs_combo=Q.QComboBox();rb.addWidget(self.obs_combo,1)
        self.mode=Q.QComboBox();self.mode.addItems(['叠加曲线','差值：本组 − 参考组','绝对差值','相对差值（除以参考值绝对值）']);rb.addWidget(self.mode)
        button('打开已有结果',self.open_results,rb);button('刷新结果',self.reload_results,rb);button('导出当前图',self.export_current,rb);button('导出全部结果',self.export_all,rb)
        self.result_context=Q.QLabel('计算后在这里显示结果');self.result_context.setTextFormat(QtCore.Qt.PlainText);rl.addWidget(self.result_context)
        self.draw_timer=QtCore.QTimer(self);self.draw_timer.setSingleShot(True);self.draw_timer.timeout.connect(self.draw)
        self.drawing=False;self.last_gui_error=None;self.last_plot=None
        self.obs_combo.currentIndexChanged.connect(self.request_draw);self.mode.currentIndexChanged.connect(self.request_draw)
        self.visible=Q.QListWidget();self.visible.setMaximumHeight(60);self.visible.itemChanged.connect(self.request_draw);rl.addWidget(self.visible)
        tools=Q.QHBoxLayout();rl.addLayout(tools)
        self.figure=Figure(figsize=(10,6));self.axes=self.figure.add_subplot(111)
        self.canvas=ResultCanvas(self.figure);self.canvas.failed.connect(self.report_gui_error)
        self.toolbar=NavigationToolbar2QT(self.canvas,self);tools.addWidget(self.toolbar,1)
        tools.addWidget(Q.QLabel('图形大小'));self.chart_zoom=Q.QSpinBox();self.chart_zoom.setRange(50,250);self.chart_zoom.setSuffix(' %');self.chart_zoom.setValue(100);tools.addWidget(self.chart_zoom)
        button('适应窗口',self.fit_chart,tools)
        self.chart_scroll=Q.QScrollArea();self.chart_scroll.setAlignment(QtCore.Qt.AlignCenter);self.chart_scroll.setWidget(self.canvas)
        self.chart_scroll.setMinimumHeight(180);self.chart_scroll.viewport().installEventFilter(self);self.canvas.installEventFilter(self)
        self.chart_zoom.valueChanged.connect(self.resize_chart)
        self.chart_resize_timer=QtCore.QTimer(self);self.chart_resize_timer.setSingleShot(True);self.chart_resize_timer.timeout.connect(self.resize_chart)
        self.result_splitter=Q.QSplitter(QtCore.Qt.Vertical);self.result_splitter.setChildrenCollapsible(False)
        self.result_splitter.addWidget(self.chart_scroll);rl.addWidget(self.result_splitter,1)
        self.result_table=Q.QTableWidget(0,6);self.result_table.setHorizontalHeaderLabels(['U','波函数组','物理量','值','状态','说明'])
        self.result_table.setMinimumHeight(95);self.result_table.setEditTriggers(Q.QAbstractItemView.NoEditTriggers)
        self.result_table.horizontalHeader().setStretchLastSection(True);self.result_splitter.addWidget(self.result_table)
        self.result_splitter.setSizes([500,150]);self.result_splitter.setStretchFactor(0,3);self.result_splitter.setStretchFactor(1,1)
        logs=Q.QWidget();self.tabs.addTab(logs,'运行记录');ll=Q.QVBoxLayout(logs)
        self.log=Q.QPlainTextEdit();self.log.setReadOnly(True);self.log.setMaximumBlockCount(5000);ll.addWidget(self.log)
        button('打开当前结果目录',self.open_run_dir,ll)
        footer=Q.QHBoxLayout();outer.addLayout(footer);self.force=Q.QCheckBox('强制重算');footer.addWidget(self.force)
        self.start_button=button('开始 / 继续计算',self.start_job,footer);self.stop_button=button('停止本次计算',self.stop_job,footer);self.stop_button.setEnabled(False)
        self.progress=Q.QProgressBar();footer.addWidget(self.progress,1);self.status=Q.QLabel('准备就绪');footer.addWidget(self.status)
        self.process=QtCore.QProcess(self);self.process.readyReadStandardOutput.connect(self.read_stdout)
        self.process.readyReadStandardError.connect(self.read_stderr);self.process.finished.connect(self.process_finished)
        self.process.errorOccurred.connect(self.process_error)
        self.save_timer=QtCore.QTimer(self);self.save_timer.setSingleShot(True);self.save_timer.timeout.connect(self.autosave)
        self.refresh_all()
        self.setStyleSheet('QMainWindow{background:#f3f6fa} QPushButton{padding:6px 10px;border:1px solid #bcc9d6;border-radius:4px;background:white} QPushButton:hover{background:#e5edf7} QTabWidget::pane{border:1px solid #c7d1dc} QHeaderView::section{background:#edf2f8;padding:5px;border:0px} QProgressBar{text-align:center} QProgressBar::chunk{background:#367ac1}')
        self.start_button.setStyleSheet('background:#2166ac;color:white;font-weight:bold;padding:8px 18px')
        self.resize_chart()
    def notify(self,e): Q.QMessageBox.warning(self,'需要处理',str(e))
    def show_definitions(self):
        try:
            if not hasattr(self,'definitions_dialog'):self.definitions_dialog=DefinitionsDialog(self)
            self.definitions_dialog.show();self.definitions_dialog.raise_();self.definitions_dialog.activateWindow()
        except Exception as e:self.notify(e)
    def report_gui_error(self,detail):
        self.status.setText('显示失败，可切换后重试；详情见运行记录')
        if detail==self.last_gui_error:return
        self.last_gui_error=detail
        try:
            with open(str(ROOT/'data'/'界面错误.log'),'a',encoding='utf-8') as f:f.write('\n'+datetime.now().isoformat()+'\n'+detail)
            self.append_log('界面显示异常（结果已保留）：\n'+detail)
        except Exception:pass
    def request_draw(self,*args):self.draw_timer.start(40)
    def fit_chart(self):
        self.chart_zoom.setValue(100);self.resize_chart()
    def resize_chart(self,*args):
        view=self.chart_scroll.viewport();scale=self.chart_zoom.value()/100.0
        width=max(640,min(view.width()-8,(view.height()-8)*1.7))
        self.canvas.setFixedSize(round(width*scale),round(width*scale/1.7))
    def eventFilter(self,obj,event):
        if hasattr(self,'chart_scroll') and obj is self.chart_scroll.viewport() and event.type()==QtCore.QEvent.Resize:
            self.chart_resize_timer.start(80)
        if hasattr(self,'canvas') and obj is self.canvas and event.type()==QtCore.QEvent.Wheel:
            horizontal=bool(event.modifiers() & QtCore.Qt.ShiftModifier)
            bar=self.chart_scroll.horizontalScrollBar() if horizontal else self.chart_scroll.verticalScrollBar()
            delta=event.pixelDelta().y() or event.angleDelta().y()/3
            bar.setValue(bar.value()-round(delta));event.accept();return True
        return super().eventFilter(obj,event)
    def autosave(self): atomic_json(ROOT/'data'/'last_project.json',self.project)
    def changed(self):self.save_timer.start(250)
    def active_group(self):
        i=self.groups.currentIndex();return self.project['groups'][max(0,i)]
    def refresh_all(self):
        active=self.groups.currentData()
        for w,v in zip(self.grid_inputs,self.project['grid']):w.setText(v)
        self.order.blockSignals(True);self.order.setValue(self.project['common']['max_order']);self.order.blockSignals(False)
        self.groups.blockSignals(True);self.groups.clear()
        self.reference.blockSignals(True);self.reference.clear()
        for g in self.project['groups']:self.groups.addItem(g['name'],g['id']);self.reference.addItem(g['name'],g['id'])
        ids=[g['id'] for g in self.project['groups']]
        if active in ids:self.groups.setCurrentIndex(ids.index(active))
        self.reference.setCurrentIndex(next(i for i,g in enumerate(self.project['groups']) if g['id']==self.project['reference']))
        self.groups.blockSignals(False);self.reference.blockSignals(False);self.refresh_tree();self.refresh_table();self.changed()
        if self.snapshot is not None and not self.running:
            current={g['id']:g for g in self.project['groups']}
            for g in self.snapshot['groups']:
                if g['id'] in current:
                    for key in ('name','color','style','notes'):g[key]=current[g['id']][key]
            self.reload_results()
    def refresh_tree(self):
        self.tree.blockSignals(True);self.tree.clear();parents={}
        for category,key,label,_ in observable_specs(self.order.value()):
            if category not in parents:
                p=Q.QTreeWidgetItem(self.tree,[category]);p.setFlags(p.flags()|QtCore.Qt.ItemIsTristate|QtCore.Qt.ItemIsUserCheckable);parents[category]=p
            item=Q.QTreeWidgetItem(parents[category],[label]);item.setData(0,QtCore.Qt.UserRole,key)
            item.setFlags(item.flags()|QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(0,QtCore.Qt.Checked if key in self.project['observables'] else QtCore.Qt.Unchecked)
        self.tree.expandAll();self.tree.blockSignals(False)
    def selection_changed(self,*args):
        keys=[]
        for i in range(self.tree.topLevelItemCount()):
            parent=self.tree.topLevelItem(i)
            for j in range(parent.childCount()):
                item=parent.child(j)
                if item.checkState(0)==QtCore.Qt.Checked:keys.append(item.data(0,QtCore.Qt.UserRole))
        self.project['observables']=keys;self.changed()
    def change_order(self,value):
        self.project['common']['max_order']=value
        valid={s[1] for s in observable_specs(value)};self.project['observables']=[k for k in self.project['observables'] if k in valid]
        self.refresh_tree();self.changed()
    def select_preset(self,kind):
        specs=observable_specs(self.order.value())
        keys=[s[1] for s in specs] if kind=='all' else [] if kind=='none' else ['half_chain_entropy','entanglement_gap','string_order','kappa2','fidelity'] if kind=='quick' else [s[1] for s in specs if s[0] in ('纠缠量','非局域量','普通累积量','标准化累积量')]
        self.project['observables']=keys;self.refresh_tree();self.changed()
    def refresh_table(self,*args):
        if not self.project['groups']:return
        g=self.active_group();us=group_grid(self.project,g);self.table.setRowCount(len(us))
        for i,u in enumerate(us):
            path=g['files'].get(u,'');status=self.checked_status.get((g['id'],u),'未检查' if path else '空缺')
            for j,value in enumerate([u,Path(path).name if path else '← 拖入波函数',path,status]):
                item=Q.QTableWidgetItem(value);item.setToolTip(value);self.table.setItem(i,j,item)
        self.pending.clear();self.pending.addItems(g['pending'])
        p=effective_params(self.project,g)
        self.group_note.setText('%s｜t1=%g，ratio=%g，dimer=%s，nbar=%g｜%s'%(g['name'],p['t1'],p['ratio'],p['dimer'],p['nbar'],'独立 U 网格' if g['grid'] else '公共 U 网格'))
    def apply_global_grid(self):
        try:update_grid(self.project,[w.text().strip() for w in self.grid_inputs]);self.checked_status.clear();self.refresh_table();self.changed()
        except Exception as e:self.notify(e)
    def import_files(self,paths,row=-1):
        try:
            g=self.active_group();us=group_grid(self.project,g)
            if len(paths)==1 and row>=0 and g['files'].get(us[row]):replace_file(g,us[row],paths[0])
            else:assign_files(self.project,g,paths,us[row] if row>=0 else None)
            self.checked_status.clear();self.refresh_table();self.changed()
        except Exception as e:self.notify(e)
    def choose_files(self):
        paths,_=Q.QFileDialog.getOpenFileNames(self,'按希望的顺序选择；导入后可交换位置','','波函数 (*.h5 *.hdf5)')
        self.import_files(paths,self.table.currentRow())
    def batch_import(self):
        try:
            g=self.active_group();us=group_grid(self.project,g);row=self.table.currentRow()
            d=BatchImportDialog(self.project,g,self,us[row] if row>=0 else None)
            if d.exec_()==Q.QDialog.Accepted:
                assign_files(self.project,g,d.paths(),d.start.currentText())
                self.checked_status.clear();self.refresh_table();self.changed()
        except Exception as e:self.notify(e)
    def replace_selected(self):
        row=self.table.currentRow()
        if row<0:return self.notify('先选择要替换的 U 行')
        path,_=Q.QFileDialog.getOpenFileName(self,'替换文件','','波函数 (*.h5 *.hdf5)')
        if path:self.import_files([path],row)
    def selected_rows(self):return sorted({i.row() for i in self.table.selectionModel().selectedRows()})
    def clear_selected(self):
        g=self.active_group();us=group_grid(self.project,g)
        for row in self.selected_rows():
            p=g['files'].pop(us[row],None)
            if p:g['pending'].append(p)
        self.checked_status.clear();self.refresh_table();self.changed()
    def swap_rows(self,a,b):
        g=self.active_group();us=group_grid(self.project,g);x=g['files'].pop(us[a],None);y=g['files'].pop(us[b],None)
        if x:g['files'][us[b]]=x
        if y:g['files'][us[a]]=y
        self.checked_status.clear();self.refresh_table();self.changed()
    def swap_selected(self):
        rows=self.selected_rows()
        if len(rows)!=2:return self.notify('请选择两行（按住 Ctrl 可多选）')
        self.swap_rows(*rows)
    def assign_pending(self):
        items=self.pending.selectedItems();paths=[x.text() for x in items]
        if not paths:return
        try:
            g=self.active_group();row=self.table.currentRow();us=group_grid(self.project,g)
            assign_files(self.project,g,paths,us[row] if row>=0 else None)
            for p in paths:g['pending'].remove(p)
            self.refresh_table();self.changed()
        except Exception as e:self.notify(e)
    def remove_pending(self):
        g=self.active_group()
        for item in self.pending.selectedItems():g['pending'].remove(item.text())
        self.refresh_table();self.changed()
    def add_group(self):self.project['groups'].append(new_group(len(self.project['groups'])));self.refresh_all();self.groups.setCurrentIndex(len(self.project['groups'])-1)
    def remove_group(self):
        if len(self.project['groups'])==1:return self.notify('至少保留一组')
        i=self.groups.currentIndex();g=self.project['groups'].pop(i)
        if g['id']==self.project['reference']:self.project['reference']=self.project['groups'][0]['id']
        self.refresh_all()
    def group_dialog(self):
        g=self.active_group();d=ParamsDialog(effective_params(self.project,g),self,overrides=bool(g['params']) or None,grid_values=g['grid'])
        form=d.inputs['dataset'].parentWidget().layout()
        namew=line(g['name']);notew=line(g['notes']);stylew=Q.QComboBox();stylew.addItems(list(STYLE));stylew.setCurrentIndex(list(STYLE.values()).index(g['style']))
        colorw=line(g['color']);colorb=Q.QPushButton('选颜色')
        def choose_color():
            c=Q.QColorDialog.getColor(QtGui.QColor(colorw.text()),self)
            if c.isValid():colorw.setText(c.name())
        colorb.clicked.connect(choose_color);cl=Q.QHBoxLayout();cl.addWidget(colorw);cl.addWidget(colorb)
        form.insertRow(0,'显示名称',namew);form.insertRow(1,'计算参数备注',notew);form.insertRow(2,'线型',stylew);form.insertRow(3,'颜色',cl)
        if d.exec_()==Q.QDialog.Accepted:
            if not namew.text().strip() or not QtGui.QColor(colorw.text()).isValid():return self.notify('名称或颜色无效')
            g['name']=namew.text().strip();g['notes']=notew.text();g['style']=STYLE[stylew.currentText()];g['color']=QtGui.QColor(colorw.text()).name()
            g['params']=d.values() if d.override.isChecked() else {}
            if d.custom.isChecked():update_grid(self.project,[w.text().strip() for w in d.ginputs],g)
            elif g['grid']:
                update_grid(self.project,self.project['grid'],g);g['grid']=None
            self.checked_status.clear();self.refresh_all()
    def common_params(self):
        d=ParamsDialog(self.project['common'],self)
        if d.exec_()==Q.QDialog.Accepted:self.project['common'].update(d.values());self.checked_status.clear();self.refresh_table();self.changed()
    def change_reference(self,index):
        if index>=0:self.project['reference']=self.reference.itemData(index);self.changed()
    def new(self):self.project=new_project();self.checked_status.clear();self.refresh_all()
    def load_example(self):self.project=example_project();self.checked_status.clear();self.refresh_all()
    def save_project(self):
        path,_=Q.QFileDialog.getSaveFileName(self,'保存对比方案',str(ROOT/'data'/'对比方案.json'),'对比方案 (*.json)')
        if path:
            if Path(path).suffix.lower()!='.json':path+='.json'
            try:atomic_json(path,self.project);self.status.setText('方案已保存')
            except Exception as e:self.notify(e)
    def open_project(self):
        path,_=Q.QFileDialog.getOpenFileName(self,'打开对比方案',str(ROOT/'data'),'对比方案 (*.json)')
        if path:
            try:self.project=validate_project(json.loads(Path(path).read_text(encoding='utf-8-sig')));self.checked_status.clear();self.refresh_all()
            except Exception as e:self.notify(e)
    def runtime_dialog(self):
        d=Q.QDialog(self);d.setWindowTitle('运行环境');layout=Q.QFormLayout(d);inputs={}
        for key,label in [('python','Python 可执行文件（下次启动生效）'),('julia','Julia 可执行文件'),('julia_project','Julia 项目目录（空 = 现有默认环境）')]:
            inputs[key]=line(self.runtime[key]);layout.addRow(label,inputs[key])
        bb=Q.QDialogButtonBox(Q.QDialogButtonBox.Ok|Q.QDialogButtonBox.Cancel);layout.addRow(bb);bb.accepted.connect(d.accept);bb.rejected.connect(d.reject)
        if d.exec_()==Q.QDialog.Accepted:
            values={k:w.text().strip() for k,w in inputs.items()}
            if not all(os.path.isfile(values[k]) for k in ('python','julia')):return self.notify('Python 或 Julia 可执行文件不存在')
            self.runtime=values;atomic_json(ROOT/'data'/'runtime.json',values)
    def check_mapping(self,notify=True):
        problems=[];fatal=False
        for g in self.project['groups']:
            p=effective_params(self.project,g)
            for u,path in g['files'].items():
                try:self.checked_status[(g['id'],u)]=check_file(path,u,p)
                except Exception as e:
                    self.checked_status[(g['id'],u)]=str(e);problems.append('%s / U=%s：%s'%(g['name'],u,e))
                    fatal=fatal or isinstance(e,MappingConflict)
        self.refresh_table()
        if problems and (notify or fatal):self.notify('\n'.join(problems))
        elif problems:self.append_log('下列文件将记录为失败，其他文件继续：\n'+'\n'.join(problems))
        else:self.status.setText('已导入文件的内部参数检查通过；空缺点将留空')
        return not fatal
    def start_job(self):
        if self.running:return
        try:
            update_grid(self.project,[w.text().strip() for w in self.grid_inputs]);self.refresh_table()
            if not self.check_mapping(False):return
            stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f');self.run_dir=ROOT/'data'/'runs'/stamp
            self.snapshot=copy.deepcopy(self.project)
            self.jobfile,points=write_job(self.snapshot,self.run_dir,self.force.isChecked())
            self.rows=[];self.completed=0;self.total=len(points)*(2 if 'center_entropy' in self.snapshot['observables'] else 1)
            self.progress.setRange(0,self.total);self.progress.setValue(0)
            self.logfile=open(str(self.run_dir/'运行日志.txt'),'a',encoding='utf-8')
            self.running=True;self.start_button.setEnabled(False);self.stop_button.setEnabled(True);self.setup.setEnabled(False)
            self.phase='normal';self.launch_worker()
        except Exception as e:self.notify(e)
    def launch_worker(self):
        args=['--startup-file=no']
        if self.runtime['julia_project']:args.append('--project='+self.runtime['julia_project'])
        args += [str(ROOT/'backend.jl'),self.jobfile,self.phase]
        self.buffer='';self.status.setText('正在启动 Julia；首次运行需要编译');self.process.setWorkingDirectory(str(ROOT))
        self.process.start(self.runtime['julia'],args)
    def append_log(self,text):
        self.log.appendPlainText(text.rstrip())
        if self.logfile:self.logfile.write(text+'\n');self.logfile.flush()
    def read_stdout(self):
        self.buffer+=bytes(self.process.readAllStandardOutput()).decode('utf-8','replace')
        while '\n' in self.buffer:
            row,self.buffer=self.buffer.split('\n',1);self.handle_event(row.rstrip('\r'))
    def handle_event(self,row):
        if not row.startswith('WC\t'):self.append_log(row);return
        fields=row.split('\t',5)
        if len(fields)!=6:self.append_log(row);return
        _,kind,g,u,key,encoded=fields
        try:message=base64.b64decode(encoded).decode('utf-8')
        except Exception:message=encoded
        names={x['id']:x['name'] for x in self.snapshot['groups']};label=names.get(g,g)
        quantity={s[1]:s[2] for s in observable_specs(32)}.get(key,key)
        human='%s｜U=%s｜%s｜%s'%(label,u,quantity,STATUS.get(kind,kind))
        if message:human+='：'+message
        self.append_log(human)
        if kind=='point_done':self.completed+=1;self.progress.setValue(self.completed);self.reload_results()
        elif kind in ('point','calculating'):self.status.setText('%s / U=%s / %s'%(label,u,quantity or ('中心块熵' if message=='center' else '检查波函数')))
        elif kind in STATUS:self.reload_results()
    def read_stderr(self):self.append_log(bytes(self.process.readAllStandardError()).decode('utf-8','replace'))
    def process_finished(self,code,status):
        self.read_stdout();self.read_stderr()
        if self.running and code==0 and self.phase=='normal' and 'center_entropy' in self.snapshot['observables']:
            self.phase='center';self.launch_worker();return
        was_running=self.running;self.running=False;self.start_button.setEnabled(True);self.stop_button.setEnabled(False);self.setup.setEnabled(True)
        self.reload_results();self.status.setText('计算结束' if was_running and code==0 else '已停止 / 可继续计算未完成项')
        if self.logfile:self.logfile.close();self.logfile=None
        if was_running and code!=0:self.notify('计算进程异常退出；已完成结果保留。详情见运行记录。')
        self.tabs.setCurrentIndex(1)
    def process_error(self,error):
        self.append_log('进程错误：'+self.process.errorString())
        if error==QtCore.QProcess.FailedToStart:
            self.running=False;self.start_button.setEnabled(True);self.stop_button.setEnabled(False);self.setup.setEnabled(True)
            if self.logfile:self.logfile.close();self.logfile=None
            self.notify('无法启动 Julia，请在运行环境中检查路径。')
    def stop_job(self):
        self.running=False;self.status.setText('正在停止本软件的计算进程');self.process.kill()
    def reload_results(self):
        if not self.run_dir or self.snapshot is None:return
        try:self.rows=read_results(self.run_dir/'results.csv')
        except (PermissionError,OSError):return
        names={g['id']:g['name'] for g in self.snapshot['groups']}
        self.result_context.setText('本次结果参考组：%s · 数值与物理设置来自开始计算时的方案快照'%names[self.snapshot['reference']])
        available=set(r['observable'] for r in self.rows)
        keys=[k for k in self.snapshot['observables'] if k in available];old=self.obs_combo.currentData()
        self.obs_combo.blockSignals(True);self.obs_combo.clear();spec={s[1]:s[2] for s in observable_specs(32)}
        for k in keys:self.obs_combo.addItem(spec.get(k,k),k)
        if old in keys:self.obs_combo.setCurrentIndex(keys.index(old))
        self.obs_combo.blockSignals(False)
        if self.visible.count()==0 or {self.visible.item(i).data(QtCore.Qt.UserRole) for i in range(self.visible.count())}!={g['id'] for g in self.snapshot['groups']}:
            self.visible.blockSignals(True);self.visible.clear()
            for g in self.snapshot['groups']:
                item=Q.QListWidgetItem(g['name']);item.setData(QtCore.Qt.UserRole,g['id']);item.setFlags(item.flags()|QtCore.Qt.ItemIsUserCheckable);item.setCheckState(QtCore.Qt.Checked);self.visible.addItem(item)
            self.visible.blockSignals(False)
        else:
            self.visible.blockSignals(True)
            names={g['id']:g['name'] for g in self.snapshot['groups']}
            for i in range(self.visible.count()):self.visible.item(i).setText(names[self.visible.item(i).data(QtCore.Qt.UserRole)])
            self.visible.blockSignals(False)
        self.draw()
    def open_results(self):
        if self.running:return self.notify('当前正在计算，请先停止或等待完成')
        directory=Q.QFileDialog.getExistingDirectory(self,'打开包含 results.csv 与方案快照的结果目录',str(ROOT/'data'/'runs'))
        if directory:
            try:
                root=Path(directory)
                self.snapshot=validate_project(json.loads((root/'方案快照.json').read_text(encoding='utf-8')))
                if not (root/'results.csv').exists():raise ValueError('没有 results.csv')
                self.project=copy.deepcopy(self.snapshot);self.run_dir=root;self.refresh_all();self.reload_results()
            except Exception as e:self.notify(e)
    def result_value(self,g,u,k):
        matches=[r for r in self.rows if r['group_id']==g and r['U']==u and r['observable']==k]
        if not matches:return float('nan')
        r=matches[-1]
        return float(r['value']) if r['status'] in ('computed','cached') else float('nan')
    def plot_axis(self,ax,key,mode=0,visible=None):
        project=self.snapshot;spec={s[1]:(s[2],s[3]) for s in observable_specs(32)}
        for g in project['groups']:
            if visible is not None and g['id'] not in visible:continue
            us=group_grid(project,g);ys=[]
            for u in us:
                value=self.result_value(g['id'],u,key)
                if mode:
                    reference=self.result_value(project['reference'],u,key)
                    value=value-reference
                    if mode==2:value=abs(value)
                    if mode==3:value=value/abs(reference) if abs(reference)>1e-12 else float('nan')
                ys.append(value)
            ax.plot([float(u) for u in us],ys,g['style'],marker='o',markersize=4,linewidth=1.6,color=g['color'],label=g['name'])
        title,ylabel=spec.get(key,(key,key));ax.set_title(title,fontsize=11);ax.set_xlabel('U');ax.set_ylabel(ylabel if not mode else ['','delta '+ylabel,'|delta '+ylabel+'|','relative delta '+ylabel][mode])
        ax.grid(alpha=.22);ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useOffset=False)
        if mode:ax.axhline(0,color='#8794a5',linewidth=.7)
        ticks=sorted({float(u) for g in project['groups'] for u in group_grid(project,g)})
        if len(ticks)<=15:ax.set_xticks(ticks)
        if key=='total_number' and mode==0:
            values=[self.result_value(g['id'],u,key) for g in project['groups'] for u in group_grid(project,g)]
            values=[v for v in values if math.isfinite(v)]
            if values and max(values)-min(values)<1e-8:
                center=sum(values)/len(values);ax.set_ylim(center-.01,center+.01)
        if ax.lines and any(l.get_label() and not l.get_label().startswith('_') for l in ax.lines):ax.legend(fontsize=8,loc='best')
    def draw(self,*args):
        if not hasattr(self,'figure') or self.snapshot is None:return
        if self.drawing:self.request_draw();return
        self.draw_timer.stop();self.drawing=True
        key=self.obs_combo.currentData()
        try:
            visible={self.visible.item(i).data(QtCore.Qt.UserRole) for i in range(self.visible.count()) if self.visible.item(i).checkState()==QtCore.Qt.Checked}
            # Keep the axes alive: the navigation toolbar holds references to them.
            identity=(key,self.mode.currentIndex(),self.run_dir)
            if identity!=self.last_plot:
                if self.toolbar.mode:
                    if 'pan' in str(self.toolbar.mode).lower():self.toolbar.pan()
                    else:self.toolbar.zoom()
                self.toolbar.update();self.last_plot=identity
            self.axes.clear();self.result_table.clearContents()
            if key:
                self.plot_axis(self.axes,key,self.mode.currentIndex(),visible)
                rows=[r for r in self.rows if r['observable']==key];names={g['id']:g['name'] for g in self.snapshot['groups']}
                self.result_table.setRowCount(len(rows))
                for i,r in enumerate(rows):
                    for j,v in enumerate([r['U'],names.get(r['group_id'],r['group_id']),self.obs_combo.currentText(),r['value'],STATUS.get(r['status'],r['status']),r['message']]):
                        item=Q.QTableWidgetItem(v);item.setToolTip(v);self.result_table.setItem(i,j,item)
            else:self.result_table.setRowCount(0)
            self.figure.subplots_adjust(left=.15,right=.96,bottom=.15,top=.9)
            self.canvas.draw_idle()
            if self.status.text().startswith('显示失败'):self.status.setText('图形显示已恢复')
        except Exception:self.report_gui_error(traceback.format_exc())
        finally:self.drawing=False
    def export_current(self):
        if not self.rows:return self.notify('先计算或打开已有结果')
        self.draw()
        path,kind=Q.QFileDialog.getSaveFileName(self,'导出当前图',str(self.run_dir/(safe_name(self.obs_combo.currentText())+'.png')),'PNG (*.png);;SVG (*.svg)')
        if path:
            if Path(path).suffix.lower() not in ('.png','.svg'):path+='.svg' if kind.startswith('SVG') else '.png'
            try:self.figure.savefig(path,dpi=220,bbox_inches='tight');self.status.setText('图片已导出')
            except Exception as e:self.notify(e)
    def export_all(self,directory=None):
        if not self.rows:return self.notify('先计算或打开已有结果')
        if not directory:directory=Q.QFileDialog.getExistingDirectory(self,'选择导出目录',str(self.run_dir))
        if not directory:return
        try:
            out=Path(directory)/('导出_'+datetime.now().strftime('%Y%m%d_%H%M%S'));out.mkdir(parents=True,exist_ok=True)
            keys=list(dict.fromkeys(r['observable'] for r in self.rows));spec={s[1]:s[2] for s in observable_specs(32)}
            for k in keys:
                fig=Figure(figsize=(9,5));ax=fig.add_subplot(111);self.plot_axis(ax,k);fig.tight_layout()
                for ext in ('png','svg'):fig.savefig(str(out/(safe_name(spec.get(k,k))+'.'+ext)),dpi=200)
            for page in range(0,len(keys),9):
                part=keys[page:page+9];fig=Figure(figsize=(15,4.3*math.ceil(len(part)/3)))
                for i,k in enumerate(part):self.plot_axis(fig.add_subplot(math.ceil(len(part)/3),3,i+1),k)
                fig.tight_layout()
                for ext in ('png','svg'):fig.savefig(str(out/('总览_%02d.'%(page//9+1)+ext)),dpi=200)
            names={g['id']:g['name'] for g in self.snapshot['groups']}
            with open(str(out/'物理量数据.csv'),'w',encoding='utf-8-sig',newline='') as f:
                fields=['group_name','group_id','U','observable','value','status','message'];w=csv.DictWriter(f,fields);w.writeheader()
                for row in self.rows:w.writerow(dict(row,group_name=names.get(row['group_id'],row['group_id'])))
            with open(str(out/'物理量差值.csv'),'w',encoding='utf-8-sig',newline='') as f:
                w=csv.writer(f);w.writerow(['group_name','reference','U','observable','value','reference_value','difference','absolute_difference','relative_difference'])
                for g in self.snapshot['groups']:
                    if g['id']==self.snapshot['reference']:continue
                    for u in group_grid(self.snapshot,g):
                        for k in keys:
                            v=self.result_value(g['id'],u,k);r=self.result_value(self.snapshot['reference'],u,k)
                            if math.isfinite(v) and math.isfinite(r):w.writerow([g['name'],names[self.snapshot['reference']],u,k,v,r,v-r,abs(v-r),abs(v-r)/abs(r) if abs(r)>1e-12 else ''])
            atomic_json(out/'对比方案.json',self.snapshot)
            for p in self.run_dir.glob('info_*.toml'):shutil.copyfile(str(p),str(out/p.name))
            definitions=(ROOT/'物理量定义.md').read_text(encoding='utf-8')
            (out/'参数与定义说明.md').write_text(definitions+'\n\n参考组：'+names[self.snapshot['reference']]+'\n实际 U、文件对应、物理设置见对比方案；实际区间与精度检查见各 info 文件。\n',encoding='utf-8')
            self.status.setText('全部结果已导出');return out
        except Exception as e:self.notify(e)
    def open_run_dir(self):
        if self.run_dir:QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(str(self.run_dir)))
    def closeEvent(self,event):
        self.autosave()
        if self.running:self.stop_job();self.process.waitForFinished(3000)
        event.accept()

def main():
    Q.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling,True)
    app=Q.QApplication(sys.argv);app.setFont(QtGui.QFont('Microsoft YaHei',9))
    project=None;last=ROOT/'data'/'last_project.json'
    if last.exists():
        try:project=validate_project(json.loads(last.read_text(encoding='utf-8')))
        except Exception:pass
    window=MainWindow(project)
    def handle_exception(kind,value,tb):window.report_gui_error(''.join(traceback.format_exception(kind,value,tb)))
    # PyQt 5 can abort on an exception escaping a signal callback.
    sys.excepthook=handle_exception
    import faulthandler
    (ROOT/'data').mkdir(parents=True,exist_ok=True)
    app.fault_log=open(str(ROOT/'data'/'崩溃诊断.log'),'a',encoding='utf-8')
    faulthandler.enable(file=app.fault_log,all_threads=True)
    screen=app.primaryScreen()
    if screen:
        area=screen.availableGeometry()
        window.resize(min(1460,max(240,area.width()-40)),min(920,max(200,area.height()-60)))
    window.show();return app.exec_()

if __name__=='__main__':
    try:sys.exit(main())
    except Exception:
        (ROOT/'data').mkdir(exist_ok=True)
        (ROOT/'data'/'startup_error.txt').write_text(traceback.format_exc(),encoding='utf-8')
        import ctypes
        ctypes.windll.user32.MessageBoxW(None,'软件启动失败，请查看 data/startup_error.txt；检查运行环境设置。','波函数对比工作台',16)
